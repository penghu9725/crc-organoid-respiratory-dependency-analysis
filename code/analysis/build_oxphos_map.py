#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build OXPHOS_SUBMODULE_MAP.tsv from Reactome GMT + curated canonical lists.
Primary module is unique per gene; secondary modules may be multiple.
mtDNA-encoded subunits are flagged (not CRISPR-targetable)."""
import urllib.request, zipfile, io, re, collections, os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("ORGANOID_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
OUTDIR = str(Path(os.environ.get("ORGANOID_METADATA_ROOT", REPO_ROOT / "data" / "external_metadata")))
os.makedirs(OUTDIR, exist_ok=True)

GMT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
req = urllib.request.Request(GMT_URL, headers={'User-Agent': 'Mozilla/5.0'})
data = urllib.request.urlopen(req, timeout=120).read()
z = zipfile.ZipFile(io.BytesIO(data))
lines = [l for l in z.read(z.namelist()[0]).decode('utf-8', errors='replace').split('\n') if l.strip()]
byid = {}
for ln in lines:
    p = ln.split('\t')
    byid[p[1]] = (p[0], p[2:])

def gset(pid):
    return set(byid.get(pid, ('', []))[1])

CI_BIOS = gset('R-HSA-6799198')
CIII_ASM = gset('R-HSA-9865881')
CIV_ASM = gset('R-HSA-9864848')
FES = gset('R-HSA-1362409')
TCA = gset('R-HSA-71403')
CRISTAE = gset('R-HSA-8949613')
ATP20 = gset('R-HSA-163210')
IMPORT = gset('R-HSA-1268020')
MTDNA_REPL = gset('R-HSA-9913635')
PYRUV = gset('R-HSA-70268')
TRANS = gset('R-HSA-5368287')
TRANS_INIT = gset('R-HSA-163282')  # POLRMT TFAM TFB2M

# ---- curated canonical lists ----
CI_STRUCT = {'NDUFA1','NDUFA2','NDUFA3','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFA10','NDUFA11','NDUFA12','NDUFA13',
             'NDUFB1','NDUFB2','NDUFB3','NDUFB4','NDUFB5','NDUFB6','NDUFB7','NDUFB8','NDUFB9','NDUFB10','NDUFB11',
             'NDUFC1','NDUFC2','NDUFS1','NDUFS2','NDUFS3','NDUFS4','NDUFS5','NDUFS6','NDUFS7','NDUFS8',
             'NDUFV1','NDUFV2','NDUFV3','NDUFAB1',
             'MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'}
CI_ASM_CUR = {'NDUFAF1','NDUFAF2','NDUFAF3','NDUFAF4','NDUFAF5','NDUFAF6','NDUFAF7','NDUFAF8',
              'TMEM126A','TMEM126B','TMEM186','TIMMDC1','DMAC1','DMAC2','FOXRED1','ACAD9','ECSIT',
              'NUBPL','LYRM2','PYURF','SFXN4'}
CII = {'SDHA','SDHB','SDHC','SDHD','SDHAF1','SDHAF2','SDHAF3','SDHAF4'}
CIII_STRUCT = {'MT-CYB','UQCRC1','UQCRC2','UQCRB','UQCRQ','UQCRH','UQCRHL','UQCR10','UQCR11','UQCRFS1','CYC1'}
CIII_ASM_CUR = {'BCS1L','UQCC1','UQCC2','UQCC3','UQCC5','UQCC6','TTC19','LYRM7','MNF1','LETM1'}
CIV_STRUCT = {'MT-CO1','MT-CO2','MT-CO3','COX4I1','COX4I2','COX5A','COX5B','COX6A1','COX6A2','COX6B1','COX6B2','COX6C',
              'COX7A1','COX7A2','COX7A2L','COX7B','COX7B2','COX7C','COX8A','COX8C','NDUFA4','NDUFA4L2'}
CIV_ASM_CUR = {'SURF1','COX10','COX11','COX14','COX15','COX16','COX17','COX18','COX19','COX20',
               'SCO1','SCO2','PET100','PET117','COA3','COA4','COA5','COA6','COA7','COA8','CMC1','CMC2','CMC4',
               'COQ10A','COQ10B','HIGD1A','HIGD2A','PNKD','RAB5IF','SMIM20','TMEM177','TMEM223','TIMM21','TACO1'}
CV_STRUCT = {'MT-ATP6','MT-ATP8','ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D','ATP5F1E',
             'ATP5MC1','ATP5MC2','ATP5MC3','ATP5PB','ATP5PD','ATP5PO','ATP5PF','ATP5MJ','ATP5MK','ATP5MG','ATP5ME','ATP5MF'}
CV_ASM_CUR = {'ATPAF1','ATPAF2','TMEM70','DMAC2L'}
FES_CUR = {'ISCU','NFS1','LYRM4','NFU1','ISCA1','ISCA2','IBA57','BOLA1','BOLA2','BOLA3','GLRX5','GLRX2',
           'FDX1','FDX2','FDXR','FXN','HSCB','HSPA9','NUBPL','ABCB7','SLC25A28','SLC25A37'}
COQ_CUR = {'PDSS1','PDSS2','COQ2','COQ3','COQ4','COQ5','COQ6','COQ7','COQ8A','COQ8B','COQ9','COQ10A','COQ10B'}
MTFAS_CUR = {'MCAT','MECR','OXSM','CBR4','HSD17B8','HTD2','ACSF3','NDUFAB1'}
TCA_CUR = {'CS','CSKMT','ACO2','IDH1','IDH2','IDH3A','IDH3B','IDH3G','OGDH','OGDHL','DLST','KGD4',
           'SUCLA2','SUCLG1','SUCLG2','MDH1','MDH2','FH','PC','PCK2','NNT','ACAT1','SIRT3','TRAP1',
           'PDHA1','PDHA2','PDHB','DLAT','DLD','PDHX','PDK1','PDK2','PDK3','PDK4','PDP1','PDP2','MPC1','MPC2','SLC25A1'}
MTDNA_CUR = {'POLG','POLG2','TWNK','SSBP1','RNASEH1','MGME1','LIG3','TOP3A','DNA2','EXOG','TEFM','TFAM','TFB1M','TFB2M','POLRMT',
             'MPV17','DGUOK','TK2','TYMP','RRM2B','SLC25A4','TOP1MT'}
IMPORT_CORE = {'TOMM20','TOMM22','TOMM40','TOMM5','TOMM6','TOMM7','TOMM70','TOMM40L','TOMM34',
               'TIMM8A','TIMM8B','TIMM9','TIMM10','TIMM10B','TIMM13','TIMM17A','TIMM17B','TIMM21','TIMM22','TIMM23','TIMM44','TIMM50','TIMMDC1',
               'SAMM50','MTX1','MTX2','PAM16','DNAJC19','GRPEL1','GRPEL2','HSPA9','HSPD1','HSPE1','PMPCA','PMPCB','PITRM1',
               'OXA1L','GFER','CHCHD2','CHCHD3','CHCHD4','CHCHD5','CHCHD7','CHCHD10','CMC2','CMC4','COA4','COA6','ROMO1'}
MEMB_CUR = {'IMMT','CHCHD3','CHCHD6','MICOS10','MICOS13','APOO','APOOL','TMEM11','DNAJC11','OPA1',
            'TAZ','CRLS1','PRELID1','TRIAP1','STARD7','SAMM50','MTX1','MTX2'}
DYN_CUR = {'MFN1','MFN2','DNM1L','FIS1','MFF','MIEF1','MIEF2','MSTO1'}
PROT_CUR = {'LONP1','CLPP','AFG3L2','SPG7','YME1L1','OMA1','HTRA2'}

MTDNA_ENC = {'MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP6','MT-ATP8'}

# ---- module definitions: (primary, gene set, source, evidence) ----
# evidence: A=canonical structural subunit; B=established assembly factor/Reactome; C=curated canonical list; D=auxiliary
MODULES = [
    ('CI_structural',        CI_STRUCT,           'KEGG hsa00190; Reactome R-HSA-6799198; HGNC group', 'A'),
    ('CI_assembly',          CI_ASM_CUR,          'Reactome R-HSA-6799198', 'B'),
    ('CII',                  CII,                 'KEGG hsa00190; Reactome R-HSA-71403', 'A'),
    ('CIII_structural',      CIII_STRUCT,         'KEGG hsa00190; Reactome R-HSA-9865881', 'A'),
    ('CIII_assembly',        CIII_ASM_CUR,        'Reactome R-HSA-9865881', 'B'),
    ('CIV_structural',       CIV_STRUCT,          'KEGG hsa00190; Reactome R-HSA-9864848', 'A'),
    ('CIV_assembly',         CIV_ASM_CUR,         'Reactome R-HSA-9864848', 'B'),
    ('CV_structural',        CV_STRUCT,           'KEGG hsa00190; Reactome R-HSA-163210', 'A'),
    ('CV_assembly',          CV_ASM_CUR,          'curated (HGNC-group-equivalent)', 'B'),
    ('FeS_biogenesis',       FES_CUR,             'Reactome R-HSA-1362409; curated', 'A'),
    ('CoQ_biosynthesis',     COQ_CUR,             'curated (MitoCarta-equivalent)', 'B'),
    ('mtFAS',                MTFAS_CUR,           'curated', 'B'),
    ('TCA_respiratory',      TCA_CUR,             'Reactome R-HSA-71403; R-HSA-70268', 'B'),
    ('mtDNA_maintenance_expression', MTDNA_CUR,   'Reactome R-HSA-9913635; R-HSA-163282; curated', 'B'),
    ('mito_import',          IMPORT_CORE,         'Reactome R-HSA-1268020', 'B'),
    ('mito_membrane_organization', MEMB_CUR,      'Reactome R-HSA-8949613; curated', 'B'),
    ('mito_dynamics',        DYN_CUR,             'curated', 'C'),
    ('mito_proteostasis',    PROT_CUR,            'curated', 'C'),
]

# mitoribosome / translation from Reactome translation set
MRP = {g for g in TRANS if g.startswith('MRPL') or g.startswith('MRPS') or g.startswith('MRP')}
MRP_ASM = {'GTPBP5','GTPBP6','GTPBP10','ERAL1','MALSU1','RBFA','METTL15','CHCHD1'}  # mitoribosome assembly factors
AARS = {g for g in TRANS if re.search(r'(AARS2|CARS2|DARS2|EARS2|FARS2|GARS1|HARS2|IARS2|KARS1|LARS2|MARS2|NARS2|PARS2|QARS1|RARS2|SARS2|TARS2|VARS2|WARS2|YARS2|SEPSECS)$', g)}
TRANS_OTHER = (TRANS - MRP - MRP_ASM - AARS) | {'MTFMT','TUFM','TSFM','GFM1','GFM2','MTIF2','MTIF3','MTO1','GTPBP3','TACO1','LRPPRC',
                                      'TRMT10C','TRMT1','TRMT2B','TRMT5','TRMT11','TRMT44','TRUB2','PUS1','PUS7L','NSUN3','NSUN4',
                                      'METTL8','ALKBH1','ALKBH8','RPUSD3','RPUSD4','TRMU','CDK5RAP1','ELP1','ELP3',
                                      'PNPT1','ELAC2','PRORP','HSD17B10','MTPAP','FASTKD2','FASTKD5','GRSF1','DHX30','SUPV3L1','SLIRP'}
# remove genes already claimed by other primaries from translation sets
MODULES.append(('mitoribosome', MRP | MRP_ASM, 'Reactome R-HSA-5368287 (MRPL/MRPS); curated assembly factors', 'C'))
MODULES.append(('mito_translation', TRANS_OTHER, 'Reactome R-HSA-5368287; curated', 'C'))

# secondary maps (primary -> secondary module)
SECONDARY = {
    'NDUFAB1': ['mtFAS'],
    'OXA1L': ['CI_assembly', 'CIV_assembly'],
    'NUBPL': ['FeS_biogenesis'],
    'HSCB': ['CIII_assembly'],
    'HSPA9': ['CIII_assembly', 'mito_import'],
    'COA1': ['CI_assembly'],
    'TIMMDC1': ['mito_import'],
    'TIMM21': ['mito_import'],
    'SDHA': ['TCA_respiratory'], 'SDHB': ['TCA_respiratory'], 'SDHC': ['TCA_respiratory'], 'SDHD': ['TCA_respiratory'],
    'FXN': ['CIII_assembly', 'TCA_respiratory', 'mito_import'],
    'ISCU': ['CIII_assembly', 'TCA_respiratory'],
    'LYRM4': ['CIII_assembly', 'TCA_respiratory'], 'NFS1': ['CIII_assembly', 'TCA_respiratory'],
    'ISCA1': ['TCA_respiratory'], 'ISCA2': ['TCA_respiratory'],
    'TACO1': ['mito_translation'], 'LRPPRC': ['mtDNA_maintenance_expression'],
    'MTX1': ['mito_membrane_organization'], 'MTX2': ['mito_membrane_organization'],
    'SAMM50': ['mito_membrane_organization'],
    'CHCHD3': ['mito_membrane_organization'],
    'BCS1L': ['mito_import'], 'CYC1': ['mito_import'],
    'POLRMT': ['mito_translation'],
    'OPA1': ['mito_dynamics'], 'YME1L1': ['mito_dynamics'], 'OMA1': ['mito_dynamics'],
    'PMPCA': ['mito_proteostasis'], 'PMPCB': ['mito_proteostasis'], 'PITRM1': ['mito_proteostasis'],
    'DMAC2L': ['CV_structural'],
}

# ---- assemble: primary unique, secondary multi ----
assign = {}   # gene -> (primary, [secondaries], source, evidence, note)
def add_primary(gene, primary, source, ev):
    if gene not in assign:
        assign[gene] = [primary, [], source, ev, '']
    else:
        assign[gene][1].append(primary)  # conflict -> becomes secondary

for primary, gset_, src, ev in MODULES:
    for g in gset_:
        add_primary(g, primary, src, ev)
for g, secs in SECONDARY.items():
    if g in assign:
        for s in secs:
            if s not in assign[g][1]:
                assign[g][1].append(s)

# verify per-gene primary uniqueness and no self-secondary
for g, (p, secs, src, ev, note) in assign.items():
    assert p not in secs, f"self-secondary for {g}: {p}"
    assert len(set(secs)) == len(secs), f"dup secondary for {g}"

rows = []
for g, (p, secs, src, ev, note) in sorted(assign.items()):
    note_parts = []
    if g in MTDNA_ENC:
        note_parts.append('mtDNA-encoded; NOT CRISPR-targetable (no sgRNA) - exclude from dependency scoring')
    if g in ('NDUFA6', 'NDUFAF8', 'ISCA2', 'MCAT', 'MSMO1', 'TYMS', 'DHFR'):
        note_parts.append('CANDIDATE-GENE')
    if g == 'NDUFA4':
        note_parts.append('reclassified: CIV subunit (COXFA4)')
    if g == 'NDUFAB1':
        note_parts.append('acyl carrier protein shared with mtFAS')
    if g == 'SDHA':
        note_parts.append('also TCA enzyme')
    rows.append([g, p, ';'.join(secs) if secs else '.', src, ev, '; '.join(note_parts)])

header = ['gene_symbol', 'submodule_primary', 'submodule_secondary', 'source', 'evidence_level', 'notes']
out = OUTDIR + r'\OXPHOS_SUBMODULE_MAP.tsv'
with open(out, 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join(header) + '\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')

print("total genes:", len(rows))
cnt = collections.Counter(r[1] for r in rows)
for m, c in sorted(cnt.items()):
    print(f"  {m:34s} {c}")
print("saved:", out)
