# Library-stratified re-evaluation of metabolic CRISPR dependencies in colorectal and oesophageal tumour organoids

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23258833.svg)](https://doi.org/10.5281/zenodo.23258833)

## Study summary

This project re-analyses published genome-scale CRISPR dependency data from patient-derived tumour organoids to examine colorectal–oesophageal (CRC–ESCA) metabolic contrasts by screening library. The library-adjusted consensus identified OXPHOS as the only metabolic pathway differing after multiple-testing correction. A biologically curated decomposition into 13 respiratory submodules was performed **post hoc**, after the pathway-level OXPHOS result; these modules were not prespecified.

Five post hoc respiratory modules were CRC-leaning in the library-adjusted consensus and remained CRC-leaning in MinLib. The five-module pattern was not reproduced in Yusa v1.1: none of the module effects was significant after correction, and four of the five estimates leaned towards ESCA. The analyses therefore do not establish a stable CRC-specific respiratory dependency phenotype in this organoid resource. A same-line comparison across nine matched PDOs screened in distinct experiments was inconclusive; a separate 22-model DepMap subset had unresolved independence.

The organoid-versus-2D comparison is descriptive and interpreted within platform; it does not establish directly comparable dependency strength across assay platforms. Patient tumour RNA and protein measurements provide molecular-abundance context, not functional dependency validation.

## Repository contents and scope

This public repository contains the code, derived data, figures and reproducibility materials associated with the historical Zenodo **v1.0.0** archive. The current BMC Genomics manuscript presents a later library-stratified reassessment. The v1.0.0 archive does **not** contain all subsequent library-, patient- and screening-quality sensitivity results.

- [`code/`](code/): analysis and figure scripts included in the v1.0.0 repository release.
- [`data/`](data/): retained derived outputs, versioned metadata, and Supplementary Tables S1–S11.
- [`figures/final/`](figures/final/): figures retained in the v1.0.0 archive; they are the historical release files and may not reflect the current manuscript’s revised figure layout.
- [`environment/`](environment/): R and Python environment information.
- [`docs/`](docs/): source-data, workflow, portability and script-mapping documentation.
- [`reproducibility/`](reproducibility/): v1.0.0 output mapping and numerical comparison records.

For the current BMC Genomics manuscript, the post hoc LOO analysis script and library-, patient- and screening-quality sensitivity scripts are supplied as manuscript Additional files 21 and 22. They are not included in the historical v1.0.0 archive. Manuscript-specific later result tables, including Tables S12–S17, are supplied with the manuscript rather than in this repository’s `data/derived/supplementary_tables/` directory. Use the manuscript Additional files for those later results; do not infer that they are part of Zenodo v1.0.0.

## Data sources

The archived project files document the public sources used for the v1.0.0 analyses, including the published organoid Figshare release, MSigDB 2025.1.Hs, HGNC, TCGA COAD/READ, GSE132465, CPTAC-2 Colon and DepMap 24Q2. Source descriptions and access information are in [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). Upstream Sanger preprocessing code is not redistributed; retrieve GPL-3.0 files directly from Garnett-Lab/SangerOrganoidBiobank at commit `0acbbd498f686627960deac08ee2f639af67f377`. Large third-party source datasets are not redistributed.

## Reproducing the archived release

The authoritative mapping from the v1.0.0 Figures 1–7 and Supplementary Tables S1–S11 to scripts and derived outputs is [`reproducibility/REPRODUCIBILITY_MATRIX.tsv`](reproducibility/REPRODUCIBILITY_MATRIX.tsv). Publication-facing scripts in this release are indexed in [`docs/SCRIPT_MAPPING.tsv`](docs/SCRIPT_MAPPING.tsv). Environment and portability instructions are in [`environment/README.md`](environment/README.md) and [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

The current manuscript’s later sensitivity outputs are provided separately with the submission and are not claimed to be reproducible from the v1.0.0 archive alone.

## Citation and archival versions

The manuscript-associated historical archive is release **v1.0.0**, archived at [https://doi.org/10.5281/zenodo.23258832](https://doi.org/10.5281/zenodo.23258832). The concept DOI for all repository versions is [https://doi.org/10.5281/zenodo.23258833](https://doi.org/10.5281/zenodo.23258833). The v1.0.0 archive remains unchanged; its contents should not be interpreted as including later manuscript sensitivity analyses. Citation metadata for the archived release are provided in [`CITATION.cff`](CITATION.cff).

## License

Original project code is released under the MIT License. The retained `data/external_metadata/DepMap_24Q2_Model.csv` remains under CC BY 4.0. Upstream Sanger code remains under GPL-3.0 and is not redistributed here. Other external source datasets are not redistributed unless explicitly stated. See [`LICENSE_NOTES.md`](LICENSE_NOTES.md) and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
