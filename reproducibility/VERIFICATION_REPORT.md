# Portable-script verification

The three publication-facing analysis scripts were executed in isolated output directories against the versioned inputs used for the retained analysis. Retained scientific outputs were not overwritten.

| Public script | Result | Acceptance rule |
|---|---|---|
| `code/analysis/organoid_vs_2d_analysis.R` | EXACT | All three scientific TSV outputs matched exactly. |
| `code/analysis/respiratory_submodule_biomarker_analysis.R` | EXACT | All five scientific TSV outputs matched exactly. |
| `code/analysis/cptac_proteogenomic_analysis.py` | NUMERICALLY_EQUIVALENT | Absolute tolerance `1e-10`; relative tolerance `1e-8`. |

For CPTAC, two tables contained floating-point last-bit differences in raw P values and BH-FDR values. The maximum absolute difference was `1.1102230246251565e-16`; the maximum relative difference was `2.0021000048493342e-15`. All schemas, identifiers, categorical fields, classifications and conclusions matched.

Overall status: **FULLY VERIFIED**.

