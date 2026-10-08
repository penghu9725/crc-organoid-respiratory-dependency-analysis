# Reproducibility

The organoid model is the unit of inference for CRISPR dependency and organoid-internal biomarker analyses. Patient is the inferential unit for TCGA, single-cell pseudobulk and CPTAC analyses; individual cells are not treated as independent patient replicates.

Large source datasets are externally hosted and must be obtained from the versioned sources in `DATA_SOURCES.md`. In particular, the 2D comparison requires DepMap 24Q2 and its exact `Model.csv`, not the latest DepMap release. Derived publication tables and small mapping files are versioned in this repository.

The respiratory-submodule script uses random seed `20260901`; the CPTAC script uses NumPy generator seed `2404`. Bootstrap resampling counts are retained unchanged in the public scripts (2,000 for respiratory-submodule confidence intervals and 10,000 for CPTAC paired-protein confidence intervals).

Portable copies were tested against retained scientific outputs. Exact equality was required for the organoid-versus-2D and respiratory-submodule analyses. CPTAC outputs were accepted as numerically equivalent only when absolute differences were at most `1e-10` or relative differences were at most `1e-8`; observed differences were limited to floating-point last-bit variation. No classification or conclusion changed.

Set `ORGANOID_PROJECT_ROOT` to the repository root when needed. `ORGANOID_SOURCE_ROOT`, `ORGANOID_DERIVED_ROOT`, `ORGANOID_METADATA_ROOT`, `ORGANOID_FIGURE_ROOT` and `ORGANOID_LOG_ROOT` can override the corresponding locations for isolated runs.

