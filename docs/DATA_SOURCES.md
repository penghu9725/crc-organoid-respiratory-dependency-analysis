# Data sources

Large external datasets are not bundled. The scripts accept `ORGANOID_SOURCE_ROOT`; their default is `data/source/` in the repository. Small version-specific mapping files are retained in `data/external_metadata/`.

| Dataset | Purpose | DOI/accession | Release/version | Expected input or retrieval instruction |
|---|---|---|---|---|
| Organoid biobank processed release | Cohort metadata, CRISPR, RNA, mutation, CNV and covariates | [10.6084/m9.figshare.28339340](https://doi.org/10.6084/m9.figshare.28339340) | Publication-associated release | Download the Figshare files and preserve their published filenames under `data/source/00_source/figshare_28339340/`. Supplementary Tables 1 and 2 provide the cohort metadata source. |
| Garnett-Lab/SangerOrganoidBiobank | Upstream reference processing | [GitHub repository](https://github.com/Garnett-Lab/SangerOrganoidBiobank) | commit `0acbbd498f686627960deac08ee2f639af67f377` | Checkout the exact commit directly from upstream. GPL-3.0 upstream scripts are not redistributed in this MIT repository. |
| MSigDB | Metabolic gene-set definitions | MSigDB | `2025.1.Hs` | Obtain the human release from MSigDB under its terms. The mapped, derived universe is retained. |
| HGNC complete gene set | Approved symbols and Ensembl mapping | HGNC complete set | Retrieval used for the publication build | Download the HGNC complete set from HGNC. The mapped, derived universe is retained. |
| TCGA COAD/READ via UCSC Xena | Tumour/normal RNA and clinical analysis | TCGA-COAD; TCGA-READ | Xena hub `2016-01-28`; GDC survival `v41.0` | Retrieve the documented COAD/READ expression, phenotype and survival files from UCSC Xena. |
| GSE132465 | Colorectal single-cell patient context | GEO `GSE132465` | Submitter-processed release | Download the processed matrices and published annotations from GEO. |
| CPTAC-2 Colon | Paired protein and tumour RNA analysis | Zenodo `8394329` | PayneLab `cptac` 1.5.14 data release | Place the five named `.gz` matrices used by `cptac_proteogenomic_analysis.py` under `data/source/00_source/P2_4B/`. The script records their expected filenames. |
| DepMap 24Q2 dependency matrices | Matched 2D cell-line comparison | DepMap 24Q2 | 24Q2 | Obtain the 24Q2 dependency matrices; do not substitute a current release. Preserve the expected RData filenames under the Figshare source directory used by the analysis. |
| DepMap 24Q2 model annotation | Cell-line lineage annotation | [10.25452/figshare.plus.25880521.v1](https://doi.org/10.25452/figshare.plus.25880521.v1) | 24Q2; file ID `46489732` | `Model.csv` is retained at `data/external_metadata/DepMap_24Q2_Model.csv` for version-specific reproducibility. MD5: `e5b60d1ce7636542d93f45494019eded`; original licence: CC BY 4.0. |

Checksums in source-provider records should be verified after download. The retained DepMap annotation is release-specific; using another release can change lineage membership.

