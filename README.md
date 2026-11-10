# UKB Dynamic Structure Function Coupling

Reproducible analysis code for the manuscript on static and dynamic structure–function coupling in 34,067 UK Biobank participants.

This release is organized in manuscript order. Each numbered script performs one analysis stage, reads only documented inputs, and writes inspectable outputs. Participant-level UK Biobank data are not included and must remain inside an approved UK Biobank environment.

## Analysis sequence

1. `R/00_extract_ukb_fields.R` optionally extracts required fields from an authorized legacy UKB bulk export.
2. `matlab/01_compute_primary_coupling.m` calculates parcel-level SFC and DSFC.
3. `matlab/02_prepare_analysis_tables.m` creates Glasser-360 and Yeo-7 tables.
4. `matlab/03_descriptive_and_age_associations.m` generates descriptive and age-association results.
   `R/03_age_gamm_visualization.R` exports the adjusted network age curves used for visualization.
5. `matlab/04_phenotype_associations.m` tests cognitive, mental-health, and physical-health phenotypes.
6. `matlab/05_lifestyle_associations.m` tests environmental and lifestyle variables.
7. `matlab/06_genetic_associations.m` tests PRS and APOE associations.
8. `matlab/07_prepare_gwas_inputs.m` defines the unrelated imaging-genetics sample and writes PLINK input files. `shell/07_run_gwas.sh` runs the 14 GWASs, `shell/07_merge_gwas_results.sh` merges chromosomes, and `R/07_gwas_postprocess.R` produces the GWAS tables, Manhattan plots, and QQ plots.
9. `R/08_lifestyle_indirect_effects.R` fits primary and SES-adjusted smoking/income indirect-effect models.
10. `R/09_age_adprs_structural_equation_models.R` fits the age and AD-PRS structural equation models.
11. `matlab/10_robustness_and_sensitivity.m` recalculates alternative coupling definitions and BOLD properties.
    `matlab/10b_intrinsic_bold_adjustment.m` repeats age, phenotype, and AD-PRS models with matched AC1 or BOLD-variability covariates.
12. `R/10_socioeconomic_association_sensitivity.R` repeats smoking/income network associations with SES covariates.
13. `R/11_formal_sfc_dsfc_age_comparison.R` formally compares SFC and DSFC age slopes.

See [the manuscript-to-code map](docs/ANALYSIS_MAP.md) for outputs corresponding to each manuscript section.

## Primary definitions

- Windowed FC is the Pearson correlation between regional BOLD time series.
- Regional window-wise SFC is the Spearman correlation between a region's SC and FC profiles, excluding its self-connection.
- The primary whole-profile calculation retains zero-valued SC entries and negative FC values.
- SFC is the temporal mean of window-wise regional SFC.
- DSFC is the temporal standard deviation of window-wise regional SFC.
- The primary window contains 30 TRs and advances by 5 TRs.

## Setup

Clone the repository and create private configuration files:

```bash
cp config/config.example.m config/config.m
cp config/config.example.R config/config.R
```

Edit only the local paths and run settings in those two files. The populated files and all derived participant-level outputs are ignored by Git.

Required software:

- MATLAB R2025a with Statistics and Machine Learning Toolbox and Parallel Computing Toolbox
- R 4.0 with `data.table`, `ggplot2`, `ggrepel`, `patchwork`, `lavaan`, `lme4`, `lmerTest`, `emmeans`, `dplyr`, `tidyr`, and `mgcv`
- PLINK 2.0 for GWAS

The required schemas are documented in [the data contract](docs/DATA_CONTRACT.md).
The committed metadata files fix the tested families at 75 phenotypes, 56 environmental/lifestyle measures, and 90 PRSs, matching the revised supplementary tables.

## Run the primary workflow

From the repository root, run:

```matlab
run('matlab/run_primary_pipeline.m')
```

This executes stages 1–7 in manuscript order. The R stages are then run separately so their package versions and bootstrap seeds are recorded explicitly:

```bash
Rscript R/08_lifestyle_indirect_effects.R
Rscript R/09_age_adprs_structural_equation_models.R
Rscript R/10_socioeconomic_association_sensitivity.R
Rscript R/11_formal_sfc_dsfc_age_comparison.R
Rscript R/99_record_session_info.R
```

After the MATLAB table-building stage, generate the Figure 2 age curves with `Rscript R/03_age_gamm_visualization.R`.

Run sensitivity calculations with:

```matlab
run('matlab/10_robustness_and_sensitivity.m')
```

For cluster execution, adapt `slurm/submit_coupling.slurm` to the local scheduler. For GWAS, first run `matlab/07_prepare_gwas_inputs.m`, then run:

```bash
PLINK2_BIN=/path/to/plink2 \
BFILE_PATTERN=/path/to/ukb_imp_chr%d \
PHENO_FILE=derived/gwas/coupling_phenotypes.tsv \
COVAR_FILE=derived/gwas/gwas_covariates.tsv \
KEEP_FILE=derived/gwas/unrelated_imaging_keep.tsv \
OUTPUT_DIR=/path/to/chromosome/results \
bash shell/07_run_gwas.sh
```

The GWAS script restricts each chromosome to the unrelated imaging sample before applying MAF, missingness, and HWE filters. After all chromosomes finish, merge and summarize the results with:

```bash
GWAS_CHR_DIR=/path/to/chromosome/results \
GWAS_MERGED_DIR=/path/to/merged/results \
bash shell/07_merge_gwas_results.sh

GWAS_DIR=/path/to/merged/results \
GWAS_REPORT_DIR=/path/to/report \
Rscript R/07_gwas_postprocess.R
```

The PRS network screen uses Bonferroni correction across 90 PRSs × 7 networks, separately for SFC and DSFC. The GWAS report shows both the conventional threshold (`P < 5 × 10^-8`) and the study-wide threshold correcting across 14 coupling phenotypes (`P < 3.57 × 10^-9`).

## Outputs

All generated content is written below `derived/`:

```text
derived/
├── coupling/
├── tables/
├── results/
│   ├── 03_age/
│   ├── 04_phenotypes/
│   ├── 05_lifestyle/
│   ├── 06_genetics/
│   ├── 08_lifestyle_indirect_effects/
│   ├── 09_age_adprs_sem/
│   ├── 10_socioeconomic_sensitivity/
│   └── 11_formal_age_comparison/
├── gwas/
└── sensitivity/
```

Estimation and publication plotting are separated. Figures should be regenerated from the exported CSV tables and MAT files; plotting choices therefore cannot change the fitted models or correction families.

## Verification

The repository includes checks that require no UK Biobank data:

```bash
python3 tests/check_repository.py
Rscript tests/parse_all_R.R
```

GitHub Actions runs the same checks on every push and pull request. They verify the expected manuscript-ordered file set, R syntax, shell syntax, and absence of private absolute paths.

## Data access

UK Biobank participant data are available only to approved researchers through the UK Biobank Access Management System. This study used UK Biobank application 19542. No participant-level imaging, phenotype, lifestyle, or genetic data are distributed here.

## Citation

Please cite the final published article associated with this repository. The manuscript title and bibliographic details will be added after acceptance.
