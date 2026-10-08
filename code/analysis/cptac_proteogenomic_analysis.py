from pathlib import Path
import os
import gzip, hashlib, json, platform
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

REPO_DEFAULT=Path(__file__).resolve().parents[2]
ROOT=Path(os.environ.get('ORGANOID_PROJECT_ROOT', REPO_DEFAULT)).resolve()
SOURCE_ROOT=Path(os.environ.get('ORGANOID_SOURCE_ROOT', ROOT/'data'/'source')).resolve()
SRC=SOURCE_ROOT/'00_source'/'P2_4B'
OUT=Path(os.environ.get('ORGANOID_DERIVED_ROOT', ROOT/'data'/'derived'/'analysis_outputs')).resolve()
FIG=Path(os.environ.get('ORGANOID_FIGURE_ROOT', ROOT/'figures'/'generated')).resolve()
LOG=Path(os.environ.get('ORGANOID_LOG_ROOT', ROOT/'reproducibility'/'runtime_logs')).resolve()
METADATA_ROOT=Path(os.environ.get('ORGANOID_METADATA_ROOT', ROOT/'data'/'external_metadata')).resolve()
OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True); LOG.mkdir(parents=True,exist_ok=True)
MODS=['CI_structural','CI_assembly','CIV_structural','CIV_assembly','FeS_biogenesis']; RNG=np.random.default_rng(2404)
MAP=METADATA_ROOT/'OXPHOS_SUBMODULE_MAP.tsv'
mp=pd.read_csv(MAP,sep='\t'); mp=mp[mp.submodule_primary.isin(MODS)&~mp.gene_symbol.str.startswith('MT-',na=False)].copy()
genes={m:sorted(mp.loc[mp.submodule_primary.eq(m),'gene_symbol'].dropna().unique()) for m in MODS}

fn={'rna':'bcm-coad-transcriptomics-CRC-gene_RSEM_tumor_normal_UQ_log2(x+1)_BCM.txt.gz','pt':'bcm-coad-proteomics-COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt.gz','pn':'bcm-coad-proteomics-COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt.gz','map':'bcm-coad-mapping-gencode.v34.basic.annotation-mapping.txt.gz','clin':'mssm-all_cancers-clinical-clinical_Pan-cancer.May2022.tsv.gz'}
gm=pd.read_csv(SRC/fn['map'],sep='\t',usecols=['gene','gene_name']).drop_duplicates('gene'); gm['ens']=gm.gene.str.replace(r'\.\d+$','',regex=True); gm=gm.drop_duplicates('ens').set_index('ens')
def load_matrix(file,first_index=True):
 d=pd.read_csv(SRC/file,sep='\t',index_col=0); d.index=d.index.astype(str).str.replace(r'\.\d+$','',regex=True); d=d.groupby(level=0).median(numeric_only=True); d['gene_symbol']=d.index.map(gm.gene_name); d=d.dropna(subset=['gene_symbol']).groupby('gene_symbol').median(numeric_only=True); return d
rna=load_matrix(fn['rna']); pt=load_matrix(fn['pt']); pn=load_matrix(fn['pn'])
clin=pd.read_csv(SRC/fn['clin'],sep='\t'); clin=clin[(clin.tumor_code=='CO')&(clin.discovery_study!='No')].copy(); clin=clin.set_index('case_id')
stage_col='baseline/tumor_stage_pathological'; clin['stage_raw']=clin[stage_col].astype(str); clin['stage_num']=clin.stage_raw.str.extract(r'([1-4IV]+)',expand=False).replace({'I':1,'II':2,'III':3,'IV':4,'1':1,'2':2,'3':3,'4':4}); clin['stage_num']=pd.to_numeric(clin.stage_num,errors='coerce')

def rank_score(mat):
 ranks=mat.rank(axis=0,pct=True,method='average'); rows=[]
 for m in MODS:
  ev=[g for g in genes[m] if g in ranks.index and mat.loc[g].notna().any()]
  s=ranks.loc[ev].median(axis=0,skipna=True) if ev else pd.Series(dtype=float)
  rows += [{'sample_id':i,'module':m,'module_score':v,'n_genes':len(ev)} for i,v in s.items()]
 return pd.DataFrame(rows)
rs=rank_score(rna); pts=rank_score(pt); pns=rank_score(pn)
paired=sorted(set(pt.columns)&set(pn.columns)); matched_rp=sorted(set(rna.columns)&set(pt.columns))

def rb_paired(x):
 x=np.asarray(x,float); x=x[np.isfinite(x)&(x!=0)]
 if len(x)==0:return np.nan
 rr=stats.rankdata(abs(x)); return float(np.sum(rr*np.sign(x))/np.sum(rr))
def boot_median_ci(x,B=10000):
 x=np.asarray(x,float); x=x[np.isfinite(x)]
 if len(x)<2:return (np.nan,np.nan)
 vals=np.median(x[RNG.integers(0,len(x),(B,len(x)))],axis=1); return tuple(np.quantile(vals,[.025,.975]))

summary=[]
for m in MODS:
 x=pts[(pts.module==m)&pts.sample_id.isin(paired)].set_index('sample_id').module_score; y=pns[(pns.module==m)&pns.sample_id.isin(paired)].set_index('sample_id').module_score; ids=sorted(set(x.index)&set(y.index)); dif=x.loc[ids].values-y.loc[ids].values; ci=boot_median_ci(dif); w=stats.wilcoxon(dif,zero_method='wilcox',alternative='two-sided')
 summary.append(dict(modality='Protein',module=m,tumor_N=pt.shape[1],normal_N=pn.shape[1],paired_N=len(ids),tumor_median=float(np.median(x.loc[ids])),normal_median=float(np.median(y.loc[ids])),paired_median_difference_tumor_minus_normal=float(np.median(dif)),CI95_low=ci[0],CI95_high=ci[1],effect_size_matched_rank_biserial=rb_paired(dif),test='paired_Wilcoxon',raw_P=w.pvalue,evaluable=True,QC_warning=''))
 summary.append(dict(modality='RNA',module=m,tumor_N=rna.shape[1],normal_N=0,paired_N=0,tumor_median=float(rs[rs.module==m].module_score.median()),normal_median=np.nan,paired_median_difference_tumor_minus_normal=np.nan,CI95_low=np.nan,CI95_high=np.nan,effect_size_matched_rank_biserial=np.nan,test='not_evaluable_no_CPTAC_normal_RNA',raw_P=np.nan,evaluable=False,QC_warning='Frozen CPTAC processed RNA matrix contains tumor samples only; no paired-normal RNA comparison was fabricated.'))
summ=pd.DataFrame(summary); summ['BH_FDR']=np.nan; ix=summ.modality.eq('Protein'); summ.loc[ix,'BH_FDR']=stats.false_discovery_control(summ.loc[ix,'raw_P'].values); summ.to_csv(OUT/'P2_4B_CPTAC_MODULE_SUMMARY.tsv',sep='\t',index=False,na_rep='NA')

coverage=[]
for modality,mat in [('RNA_tumor',rna),('Protein_tumor',pt),('Protein_normal',pn)]:
 for m in MODS:
  ev=[g for g in genes[m] if g in mat.index and mat.loc[g].notna().any()]; coverage.append(dict(modality=modality,module=m,requested_genes=len(genes[m]),detected_genes=len(ev),coverage_fraction=len(ev)/len(genes[m]),missing_genes=';'.join(sorted(set(genes[m])-set(ev))),sample_N=mat.shape[1]))
pd.DataFrame(coverage).to_csv(OUT/'P2_4B_CPTAC_GENE_COVERAGE.tsv',sep='\t',index=False,na_rep='NA')

conc=[]
for m in MODS:
 for g in genes[m]:
  if g not in rna.index or g not in pt.index: continue
  ids=[i for i in matched_rp if pd.notna(rna.at[g,i]) and pd.notna(pt.at[g,i])]
  if len(ids)<10: continue
  rho,p=stats.spearmanr(rna.loc[g,ids],pt.loc[g,ids]); pdir=np.nan
  if g in pn.index:
   q=[i for i in paired if pd.notna(pt.at[g,i]) and pd.notna(pn.at[g,i])]; pdir=float(np.median(pt.loc[g,q].values-pn.loc[g,q].values)) if q else np.nan
  conc.append(dict(level='gene',feature=g,module=m,N=len(ids),spearman_rho=rho,raw_P=p,correlation_direction='positive_concordance' if rho>0 else 'inverse_concordance',protein_paired_direction='tumor_higher' if pdir>0 else ('tumor_lower' if pdir<0 else 'NA'),protein_paired_median_difference=pdir,QC_warning='RNA tumor-normal direction unavailable'))
 for_m_r=rs[(rs.module==m)&rs.sample_id.isin(matched_rp)].set_index('sample_id').module_score; for_m_p=pts[(pts.module==m)&pts.sample_id.isin(matched_rp)].set_index('sample_id').module_score; ids=sorted(set(for_m_r.index)&set(for_m_p.index)); rho,p=stats.spearmanr(for_m_r.loc[ids],for_m_p.loc[ids]); pdif=summ[(summ.modality=='Protein')&(summ.module==m)].paired_median_difference_tumor_minus_normal.iloc[0]
 conc.append(dict(level='module',feature=m,module=m,N=len(ids),spearman_rho=rho,raw_P=p,correlation_direction='positive_concordance' if rho>0 else 'inverse_concordance',protein_paired_direction='tumor_higher' if pdif>0 else 'tumor_lower',protein_paired_median_difference=pdif,QC_warning='RNA tumor-normal direction unavailable'))
con=pd.DataFrame(conc); con['BH_FDR']=np.nan
for lev in ['gene','module']:
 ix=con.level.eq(lev); con.loc[ix,'BH_FDR']=stats.false_discovery_control(con.loc[ix,'raw_P'].values)
con.to_csv(OUT/'P2_4B_RNA_PROTEIN_CONCORDANCE.tsv',sep='\t',index=False,na_rep='NA')

clinical=[]
for modality,sc in [('RNA',rs),('Protein',pts)]:
 for m in MODS:
  d=sc[sc.module.eq(m)].set_index('sample_id')[['module_score']].join(clin[['stage_num','stage_raw']],how='inner').dropna(subset=['module_score','stage_num']); ns=d.groupby('stage_num').size(); under=(len(ns)<2 or (ns<5).any())
  if len(d)>=10 and d.stage_num.nunique()>=2:
   rho,p=stats.spearmanr(d.stage_num,d.module_score); kw=stats.kruskal(*[x.module_score.values for _,x in d.groupby('stage_num')]); eps=max(0,(kw.statistic-d.stage_num.nunique()+1)/(len(d)-d.stage_num.nunique()))
   clinical.append(dict(modality=modality,module=m,association='Stage_Spearman',N=len(d),groups=';'.join(f'{int(k)}:{v}' for k,v in ns.items()),effect=rho,effect_label='spearman_rho',raw_P=p,underpowered=under,availability='available'))
   clinical.append(dict(modality=modality,module=m,association='Stage_Kruskal',N=len(d),groups=';'.join(f'{int(k)}:{v}' for k,v in ns.items()),effect=eps,effect_label='epsilon_squared',raw_P=kw.pvalue,underpowered=under,availability='available'))
  clinical.append(dict(modality=modality,module=m,association='MSI',N=0,groups='',effect=np.nan,effect_label='',raw_P=np.nan,underpowered=True,availability='not_available_in_frozen_CPTAC_clinical_table'))
cl=pd.DataFrame(clinical); cl['BH_FDR']=np.nan
for key,g in cl[cl.association!='MSI'].groupby(['modality','association']): cl.loc[g.index,'BH_FDR']=stats.false_discovery_control(g.raw_P.values)
cl.to_csv(OUT/'P2_4B_CPTAC_CLINICAL_ASSOCIATIONS.tsv',sep='\t',index=False,na_rep='NA')

# figures
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False})
colors={'Normal':'#56B4E9','Tumor':'#D55E00'}
fig,axes=plt.subplots(2,5,figsize=(12,5.5),constrained_layout=True)
for j,m in enumerate(MODS):
 x=pts[(pts.module==m)&pts.sample_id.isin(paired)].set_index('sample_id').module_score; y=pns[(pns.module==m)&pns.sample_id.isin(paired)].set_index('sample_id').module_score; ids=sorted(set(x.index)&set(y.index))
 ax=axes[0,j]
 for i in ids: ax.plot([0,1],[y[i],x[i]],color='0.75',lw=.5,alpha=.5)
 ax.scatter(np.zeros(len(ids)),y.loc[ids],s=7,color=colors['Normal'],alpha=.65); ax.scatter(np.ones(len(ids)),x.loc[ids],s=7,color=colors['Tumor'],alpha=.65); ax.set_xticks([0,1],['Normal','Tumor'],rotation=25); ax.set_title(m); ax.set_ylabel('Protein rank score' if j==0 else '')
 r=rs[(rs.module==m)&rs.sample_id.isin(matched_rp)].set_index('sample_id').module_score; p=pts[(pts.module==m)&pts.sample_id.isin(matched_rp)].set_index('sample_id').module_score; ids2=sorted(set(r.index)&set(p.index)); ax=axes[1,j]; ax.scatter(r.loc[ids2],p.loc[ids2],s=10,alpha=.65,color='#009E73'); rho=con[(con.level=='module')&(con.module==m)].spearman_rho.iloc[0]; ax.text(.04,.94,f'Spearman ρ={rho:.2f}\nn={len(ids2)}',transform=ax.transAxes,va='top'); ax.set_xlabel('Tumor RNA rank score'); ax.set_ylabel('Tumor protein rank score' if j==0 else '')
fig.suptitle('CPTAC Colon RNA and protein module programs\nAbundance programs are not CRISPR dependency',fontsize=12)
fig.savefig(FIG/'P2_F14_CPTAC_RNA_protein_modules.png',dpi=300); fig.savefig(FIG/'P2_F14_CPTAC_RNA_protein_modules.pdf'); plt.close(fig)

g=con[con.level=='gene'].copy(); order=g.groupby('module').spearman_rho.median().sort_values().index; fig,ax=plt.subplots(figsize=(8,7),constrained_layout=True)
pos={m:i for i,m in enumerate(order)}; jitter={m:np.linspace(-.28,.28,max(1,(g.module==m).sum())) for m in order}
for m in order:
 d=g[g.module==m].sort_values('spearman_rho'); ax.scatter(d.spearman_rho, np.array([pos[m]]*len(d))+jitter[m], c=np.where(d.BH_FDR<.05,'#0072B2','0.65'),s=16,alpha=.85,label=None)
mods_con=con[con.level=='module'].set_index('module'); ax.scatter(mods_con.loc[order].spearman_rho,range(len(order)),marker='D',s=65,color='#D55E00',edgecolor='black',label='Module-level'); ax.axvline(0,color='0.4',ls='--'); ax.set_yticks(range(len(order)),order); ax.set_xlabel('Tumor RNA–protein Spearman rho'); ax.set_title('CPTAC Colon RNA–protein concordance\nBlue genes: gene-level BH-FDR < 0.05'); ax.legend(frameon=False)
fig.savefig(FIG/'P2_F15_CPTAC_RNA_protein_concordance.png',dpi=300); fig.savefig(FIG/'P2_F15_CPTAC_RNA_protein_concordance.pdf'); plt.close(fig)

# source freeze
source_files=[fn['rna'],fn['pt'],fn['pn'],fn['map'],fn['clin'],'umich-coad-mapping-CRC_Prospective sample info.xlsx','cptac_index.tsv','cptac-master.zip']
rows=[]
for f in source_files:
 p=SRC/f; h=hashlib.sha256(p.read_bytes()).hexdigest(); rows.append(dict(file=f,bytes=p.stat().st_size,sha256=h,source='Zenodo 8394329 / PayneLab cptac 1.5.14' if f not in ['cptac-master.zip'] else 'PayneLab/cptac master snapshot',source_version='Zenodo record 8394329; cptac package 1.5.14',download_date='2026-09-02'))
pd.DataFrame(rows).to_csv(LOG/'P2_4B_SOURCE_FREEZE.tsv',sep='\t',index=False)
qc={'tumor_protein_N':pt.shape[1],'normal_protein_N':pn.shape[1],'paired_protein_N':len(paired),'tumor_RNA_N':rna.shape[1],'matched_tumor_RNA_protein_N':len(matched_rp),'clinical_N':len(clin),'stage_available_N':int(clin.stage_num.notna().sum()),'MSI_available_N':0,'module_gene_counts':{m:len(genes[m]) for m in MODS},'software':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__}}
(LOG/'P2_4B_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf-8'); print(json.dumps(qc,indent=2))
