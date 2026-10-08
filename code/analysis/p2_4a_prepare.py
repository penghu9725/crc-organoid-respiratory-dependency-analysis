import os, gzip, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(os.environ.get('ORGANOID_PROJECT_ROOT',Path(__file__).resolve().parents[2])).resolve()
SOURCE_ROOT=Path(os.environ.get('ORGANOID_SOURCE_ROOT',ROOT/'data'/'source')).resolve()
SRC=SOURCE_ROOT/'00_source'/'P2_4A'; TMP=Path(os.environ.get('ORGANOID_P2_4A_TEMP',ROOT/'data'/'source'/'.tmp_p2_4a')).resolve()
TMP.mkdir(exist_ok=True)
MODS=['CI_structural','CI_assembly','CIV_structural','CIV_assembly','FeS_biogenesis']
MAP=Path(os.environ.get('ORGANOID_METADATA_ROOT',ROOT/'data'/'external_metadata'))/'OXPHOS_SUBMODULE_MAP.tsv'
mp=pd.read_csv(MAP,sep='\t'); mp=mp[mp.submodule_primary.isin(MODS)&~mp.gene_symbol.str.startswith('MT-')].drop_duplicates(['gene_symbol','submodule_primary'])
genes=set(mp.gene_symbol)

def first_nonnull(s):
    x=s.dropna(); return x.iloc[0] if len(x) else np.nan

tcga_scores=[]; tcga_cov=[]; coverage=[]
for project in ['COAD','READ']:
    exp=pd.read_csv(SRC/f'TCGA_{project}_HiSeqV2.gz',sep='\t',index_col=0)
    exp=exp[~exp.index.duplicated(keep='first')]
    samples=[]
    for s in exp.columns:
        code=s[13:15] if len(s)>=15 else ''
        status='Tumor' if code in [f'{i:02d}' for i in range(1,10)] else ('Normal' if code in [f'{i:02d}' for i in range(10,20)] else 'Other')
        samples.append((s,s[:12],status))
    sm=pd.DataFrame(samples,columns=['sample_id','patient_id','status']); sm=sm[sm.status.isin(['Tumor','Normal'])].drop_duplicates(['patient_id','status'])
    exp=exp.loc[:,sm.sample_id]
    ranks=exp.rank(axis=0,method='average',pct=True)
    for mod in MODS:
        requested=mp.loc[mp.submodule_primary.eq(mod),'gene_symbol'].unique(); present=[g for g in requested if g in ranks.index]
        coverage.append(dict(dataset=f'TCGA-{project}',module=mod,requested_genes=len(requested),detected_genes=len(present),coverage_fraction=len(present)/len(requested),missing_genes=';'.join(sorted(set(requested)-set(present))),criterion='gene row present in HiSeqV2'))
        score=ranks.loc[present].median(axis=0)
        tcga_scores.append(pd.DataFrame({'sample_id':score.index,'module':mod,'module_score':score.values}))
    cli=pd.read_csv(SRC/f'TCGA_{project}_clinicalMatrix.tsv',sep='\t')
    cli['patient_id']=cli.sampleID.str[:12]
    cols=['patient_id','age_at_initial_pathologic_diagnosis','gender','pathologic_stage','CDE_ID_3226963']
    cli=cli[cols].groupby('patient_id',as_index=False).agg(first_nonnull)
    cli=cli.rename(columns={'age_at_initial_pathologic_diagnosis':'age','gender':'sex','pathologic_stage':'stage','CDE_ID_3226963':'MSI'})
    surv=pd.read_csv(SRC/f'TCGA_{project}_survival.tsv.gz',sep='\t'); surv['patient_id']=surv['_PATIENT'].fillna(surv['sample'].str[:12]); surv=surv.groupby('patient_id',as_index=False).agg({'OS.time':'max','OS':'max'})
    cv=sm.merge(cli,on='patient_id',how='left').merge(surv,on='patient_id',how='left'); cv['project']=project
    tcga_cov.append(cv)
tcga_scores=pd.concat(tcga_scores,ignore_index=True); tcga_cov=pd.concat(tcga_cov,ignore_index=True)
tcga_scores=tcga_scores.merge(tcga_cov,on='sample_id',how='left')
tcga_scores.to_csv(TMP/'tcga_scores.tsv',sep='\t',index=False)

# GSE132465 streaming UMI aggregation. Keep module-gene rows and total UMI per cell.
ann=pd.read_csv(SRC/'GSE132465_cell_annotation.txt.gz',sep='\t')
umi=SRC/'GSE132465_raw_UMI_count_matrix.txt.gz'
selected={}
with gzip.open(umi,'rb') as fh:
    header=fh.readline().rstrip(b'\r\n').split(b'\t'); cell_ids=np.asarray([x.decode() for x in header[1:]])
    totals=np.zeros(len(cell_ids),dtype=np.int64)
    for line in fh:
        cut=line.find(b'\t'); gene=line[:cut].decode(); vals=np.fromstring(line[cut+1:],sep='\t',dtype=np.int32)
        if len(vals)!=len(cell_ids): raise ValueError(f'UMI row width mismatch: {gene} {len(vals)} != {len(cell_ids)}')
        totals += vals
        if gene in genes and gene not in selected: selected[gene]=vals.copy()
sel=pd.DataFrame.from_dict(selected,orient='index',columns=cell_ids)
if not np.array_equal(cell_ids,ann['Index'].to_numpy()):
    pos=pd.Series(np.arange(len(cell_ids)),index=cell_ids); ann=ann[ann.Index.isin(pos.index)].copy(); order=pos.loc[ann.Index].to_numpy(); sel=sel.iloc[:,order]; totals=totals[order]

group_defs=[]
for level,keys in [('patient_class',['Patient','Class']),('patient_class_celltype',['Patient','Class','Cell_type'])]:
    codes,uniques=pd.factorize(pd.MultiIndex.from_frame(ann[keys]))
    for code,key in enumerate(uniques):
        idx=np.where(codes==code)[0]; n=len(idx)
        if n<100: continue
        lib=totals[idx].sum()
        counts=sel.iloc[:,idx].sum(axis=1).to_numpy()
        vals=np.log1p(counts/lib*1e6) if lib>0 else np.full(len(counts),np.nan)
        meta=dict(zip(keys,key if isinstance(key,tuple) else (key,))); meta.update(level=level,n_cells=n,total_UMI=int(lib))
        for gene,val,count in zip(sel.index,vals,counts): group_defs.append({**meta,'gene':gene,'logCPM':val,'sum_UMI':int(count)})
pb=pd.DataFrame(group_defs)
pb['gene_z']=pb.groupby(['level','gene']).logCPM.transform(lambda x:(x-x.mean())/x.std(ddof=0) if x.std(ddof=0)>0 else 0)
pb=pb.merge(mp[['gene_symbol','submodule_primary']].rename(columns={'gene_symbol':'gene','submodule_primary':'module'}),on='gene')
scores=pb.groupby(['level','Patient','Class','Cell_type','n_cells','total_UMI','module'],dropna=False).agg(module_score=('gene_z','median'),mean_logCPM=('logCPM','mean'),n_detected_genes=('sum_UMI',lambda x:int((x>0).sum())),n_module_genes=('gene','nunique')).reset_index()
scores.to_csv(TMP/'sc_scores.tsv',sep='\t',index=False)

for mod in MODS:
    req=mp.loc[mp.submodule_primary.eq(mod),'gene_symbol'].unique(); present=[g for g in req if g in sel.index]
    sub=sel.loc[present] if present else pd.DataFrame()
    coverage.append(dict(dataset='GSE132465',module=mod,requested_genes=len(req),detected_genes=int((sub.sum(axis=1)>0).sum()) if len(sub) else 0,coverage_fraction=(sub.sum(axis=1)>0).mean() if len(sub) else 0,missing_genes=';'.join(sorted(set(req)-set(present))),criterion='gene detected with >=1 UMI in dataset'))
pd.DataFrame(coverage).to_csv(TMP/'gene_coverage.tsv',sep='\t',index=False)

tcga_counts = tcga_cov.groupby(['project','status']).size().reset_index(name='N')
qc={'GSE132465':{'patients':int(ann.Patient.nunique()),'samples':int(ann.Sample.nunique()),'tumor_samples':int(ann.loc[ann.Class.eq('Tumor'),'Sample'].nunique()),'normal_samples':int(ann.loc[ann.Class.eq('Normal'),'Sample'].nunique()),'total_cells':len(ann),'malignant_epithelial_cells':int(((ann.Class=='Tumor')&(ann.Cell_type=='Epithelial cells')).sum()),'normal_epithelial_cells':int(((ann.Class=='Normal')&(ann.Cell_type=='Epithelial cells')).sum()),'cells_per_patient':ann.groupby('Patient').size().to_dict()},'TCGA':tcga_counts.to_dict(orient='records')}
(TMP/'qc.json').write_text(json.dumps(qc,indent=2),encoding='utf-8')

urls={
'TCGA_COAD_HiSeqV2.gz':'https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.COAD.sampleMap%2FHiSeqV2.gz','TCGA_READ_HiSeqV2.gz':'https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.READ.sampleMap%2FHiSeqV2.gz','TCGA_COAD_clinicalMatrix.tsv':'https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.COAD.sampleMap%2FCOAD_clinicalMatrix','TCGA_READ_clinicalMatrix.tsv':'https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.READ.sampleMap%2FREAD_clinicalMatrix','TCGA_COAD_survival.tsv.gz':'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.survival.tsv.gz','TCGA_READ_survival.tsv.gz':'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.survival.tsv.gz','GSE132465_cell_annotation.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz','GSE132465_raw_UMI_count_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz'}
rows=[]
for f,u in urls.items():
    p=SRC/f; h=hashlib.sha256();
    with p.open('rb') as z:
        for b in iter(lambda:z.read(8*1024*1024),b''):h.update(b)
    rows.append({'file':f,'dataset':'GSE132465' if f.startswith('GSE') else 'TCGA-COAD/READ','source_url':u,'source_version':'GEO submitter processed 2019/2020' if f.startswith('GSE') else ('UCSC TCGA Hub 2016-01-28' if 'HiSeqV2' in f or 'clinical' in f else 'GDC Xena survival v41.0'),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
pd.DataFrame(rows).to_csv(ROOT/'07_logs'/'P2_4A_SOURCE_FREEZE.tsv',sep='\t',index=False)
print(json.dumps(qc,indent=2))
