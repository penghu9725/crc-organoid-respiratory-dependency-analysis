# Analysis workflow

1. Build audited organoid metadata, CRISPR consensus/QC and the metabolic gene universe with `build_stage_p2_0_1.py`.
2. Build the prespecified respiratory-submodule map with `build_oxphos_map.py`.
3. Run the 1,412-gene and 21-pathway metabolic dependency analysis with `run_p2_stage2_scout.R`.
4. Run `organoid_vs_2d_analysis.R` for the candidate-focused organoid-versus-DepMap comparison.
5. Run `respiratory_submodule_biomarker_analysis.R` for respiratory-submodule and focused biomarker analyses.
6. Prepare and analyse TCGA COAD/READ and GSE132465 with `p2_4a_prepare.py` and `p2_4a_analysis.R`.
7. Run `cptac_proteogenomic_analysis.py` for CPTAC paired-protein and RNA-protein analyses.
8. Render Figures 1–7 using `code/figures/fig1.py` through `fig7.py` from the retained derived tables.

The output-to-script mapping is provided in `reproducibility/REPRODUCIBILITY_MATRIX.tsv`. External RNA and protein measurements represent molecular abundance, not functional dependency.

