# Manuscript to code map

The scripts follow the order of the Methods and Results in the revised manuscript. Each stage reads an explicit input contract and writes inspectable CSV files under `derived/`.

| Order | Manuscript analysis | Script | Main outputs |
|---|---|---|---|
| 0 | UKB field extraction | `R/00_extract_ukb_fields.R` | source-field subsets and manifest |
| 1 | Static and dynamic structure function coupling | `matlab/01_compute_primary_coupling.m` | participant by parcel SFC and DSFC |
| 2 | Glasser parcel and Yeo seven-network tables | `matlab/02_prepare_analysis_tables.m` | `coupling_360.csv`, `coupling_7n.csv` |
| 3 | Spatial organization and age associations | `matlab/03_descriptive_and_age_associations.m` | network descriptives and age LMMs |
| 3b | Network age curves for Figure 2C | `R/03_age_gamm_visualization.R` | GAMM predictions and fitted models |
| 4 | Cognitive, mental-health, and physical-health associations | `matlab/04_phenotype_associations.m` | network screen and parcel follow-up |
| 5 | Environmental and lifestyle associations | `matlab/05_lifestyle_associations.m` | network screen and parcel follow-up |
| 6 | PRS and APOE associations | `matlab/06_genetic_associations.m` | network PRS and parcel AD/APOE results |
| 7 | GWAS | `matlab/07_prepare_gwas_inputs.m`, `shell/07_run_gwas.sh`, `shell/07_merge_gwas_results.sh`, `R/07_gwas_postprocess.R` | PLINK inputs, chromosome and merged results, summary tables, Manhattan plots, and QQ plots |
| 8 | Smoking/income statistical indirect effects | `R/08_lifestyle_indirect_effects.R` | primary and SES-adjusted lavaan outputs |
| 9 | Age and AD-PRS structural equation models | `R/09_age_adprs_structural_equation_models.R` | SFC and DSFC lavaan outputs |
| 10 | Methodological sensitivity analyses | `matlab/10_robustness_and_sensitivity.m` | alternative estimates and correspondence tables |
| 10a | Intrinsic BOLD property adjustment | `matlab/10b_intrinsic_bold_adjustment.m` | AC1- and variability-adjusted association tables |
| 10b | SES-adjusted network associations | `R/10_socioeconomic_association_sensitivity.R` | smoking and income LMM table |
| 11 | Formal SFC versus DSFC age-slope comparison | `R/11_formal_sfc_dsfc_age_comparison.R` | omnibus test, slopes, paired contrasts |

Plotting is intentionally separated from estimation. Publication figures should be regenerated from the exported result tables so that figure styling cannot alter any model or multiple-testing result.
