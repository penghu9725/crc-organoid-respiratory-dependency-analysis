# Respiratory Complex I/IV dependencies in colorectal cancer organoids

## Overview

This repository contains the code, retained derived data and final figures supporting a comparative analysis of metabolic CRISPR dependencies in colorectal and oesophageal cancer organoids. It resolves the broad oxidative-phosphorylation signal into prespecified respiratory submodules, compares candidate dependencies with matched 2D cell lines, and evaluates RNA and protein abundance programmes in patient datasets. Patient RNA and protein abundance is not interpreted as CRISPR dependency.

## Repository structure

- `code/`: preprocessing, analysis and Figure 1–7 rendering scripts.
- `data/`: retained derived results, Supplementary Tables S1–S11 and small versioned metadata.
- `figures/final/`: unchanged publication versions of Figures 1–7 in PDF and SVG.
- `environment/`: verified R and Python environments.
- `docs/`: source-data, workflow, portability and script-mapping documentation.
- `reproducibility/`: output mapping and independent numerical comparison records.

## Data sources

Source datasets include the publication-associated organoid Figshare release, SangerOrganoidBiobank code, MSigDB 2025.1.Hs, HGNC, TCGA COAD/READ, GSE132465, CPTAC-2 Colon and DepMap 24Q2. Large third-party datasets are not redistributed. See [docs/DATA_SOURCES.md](docs/DATA_SOURCES.md).

## Analysis workflow

1. Audit organoid metadata and construct the CRISPR consensus.
2. Build the metabolic gene universe and respiratory-submodule map.
3. Test metabolic genes and pathways in CRC versus ESCA organoids.
4. Perform respiratory-submodule decomposition.
5. Compare focused candidates with lineage-matched DepMap 24Q2 cell lines.
6. Test focused organoid RNA and genomic biomarkers.
7. Evaluate patient RNA programmes in TCGA COAD/READ.
8. Evaluate patient-level single-cell programmes in GSE132465.
9. Evaluate CPTAC paired protein abundance and RNA-protein concordance.
10. Render Figures 1–7 from retained numerical source tables.

Detailed execution order is in [docs/ANALYSIS_WORKFLOW.md](docs/ANALYSIS_WORKFLOW.md).

## Reproducing the manuscript results

The authoritative mapping from Figures 1–7 and Supplementary Tables S1–S11 to source data, analysis scripts, derived tables and rendering scripts is [reproducibility/REPRODUCIBILITY_MATRIX.tsv](reproducibility/REPRODUCIBILITY_MATRIX.tsv). The publication-facing analysis scripts are listed in [docs/SCRIPT_MAPPING.tsv](docs/SCRIPT_MAPPING.tsv). Final Supplementary Tables are under `data/derived/supplementary_tables/`.

## Software environment

Verified dependencies are recorded under [environment/](environment/README.md). Public analysis scripts use repository-relative defaults and optional environment-variable overrides described in [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

## Reproducibility verification

- Organoid-versus-2D-derived outputs reproduced exactly.
- Respiratory-submodule-derived outputs reproduced exactly.
- CPTAC-derived outputs reproduced within predefined floating-point tolerance.
- No scientific classification or conclusion differed.

## Citation

The manuscript citation and DOI will be added at publication. A versioned template is provided in `CITATION.cff`.

## License

Original analysis code is released under the MIT License. Third-party datasets remain governed by their original providers and are not covered by the repository code licence; see `LICENSE_NOTES.md`.

