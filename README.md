# UKB Dynamic Structure Function Coupling

Reproducible analysis code for the manuscript "Stable brain architecture but variable dynamics: structure–function coupling signatures associated with age, genetic Alzheimer’s risk, and cognitive performance in 34,067 adults".

## Analysis sequence

1. `R/00_extract_ukb_fields.R` optionally extracts required fields from an authorized legacy UKB bulk export.
2. `matlab/01_compute_primary_coupling.m` calculates parcel-level SFC and DSFC.
3. `matlab/02_prepare_analysis_tables.m` creates Glasser-360 and Yeo-7 tables.
4. `matlab/03_descriptive_and_age_associations.m` generates descriptive and age-association results.
   `R/03_age_gamm_visualization.R` exports the adjusted network age curves used for visualization.
5. `matlab/04_phenotype_associations.m` tests cognitive, mental-health, and physical-health phenotypes.
6. `matlab/05_lifestyle_associations.m` tests environmental and lifestyle variables.
7. `matlab/06_genetic_associations.m` tests PRS and APOE associations.
8. `matlab/07_prepare_gwas_inputs.m` and `shell/07_run_gwas.sh` prepare and run GWAS.
9. `R/08_lifestyle_indirect_effects.R` fits primary and SES-adjusted smoking/income indirect-effect models.
10. `R/09_age_adprs_structural_equation_models.R` fits the age and AD-PRS structural equation models.
11. `matlab/10_robustness_and_sensitivity.m` recalculates alternative coupling definitions and BOLD properties.
    `matlab/10b_intrinsic_bold_adjustment.m` repeats age, phenotype, and AD-PRS models with matched AC1 or BOLD-variability covariates.
12. `R/10_socioeconomic_association_sensitivity.R` repeats smoking/income network associations with SES covariates.
13. `R/11_formal_sfc_dsfc_age_comparison.R` formally compares SFC and DSFC age slopes.

## Required software

- MATLAB R2025a with Statistics and Machine Learning Toolbox and Parallel Computing Toolbox; EMIGMA_toolbox (https://github.com/MICA-MNI/ENIGMA) for visualize the brain map
- R 4.2.3 with `data.table`, `lavaan`, `lme4`, `lmerTest`, `emmeans`, `dplyr`, `tidyr`, and `mgcv`
- PLINK 2.0 for GWAS

## Data access

UK Biobank participant data are available only to approved researchers through the UK Biobank Access Management System. This study used UK Biobank application 19542. No participant-level imaging, phenotype, lifestyle, or genetic data are distributed here.
