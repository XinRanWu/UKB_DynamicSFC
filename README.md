# UKB Dynamic Structure Function Coupling

Analysis code for the manuscript **Stable brain architecture but flexible dynamics: distinct structure-function coupling signatures link aging and genetic Alzheimer's risk to cognitive decline in 34,067 adults**.

The repository is organized in manuscript order. Primary analyses are executable scripts. Small MATLAB functions are limited to reusable mixed-model utilities in `utilities/`.

## What the pipeline does

1. Extracts the required UK Biobank tabular fields from an authorized local export.
2. Matches resting-state fMRI and diffusion MRI participants by `eid`.
3. Computes regional SFC in sliding windows and summarizes each participant by:
   - SFC: temporal mean of window-level regional SFC.
   - DSFC: temporal standard deviation of window-level regional SFC.
4. Builds Glasser-360 and Yeo-7 analysis tables.
5. Runs demographic, age, phenotype, lifestyle, PRS, APOE, and AD-variant mixed models.
6. Prepares network phenotypes for PLINK 2 GWAS.
7. Runs statistical indirect-effect analyses for age/APOE, coupling, and fluid intelligence.
8. Quantifies robustness to window length and structural-connectivity definition.

Participant-level UK Biobank data are not distributed in this repository. Reproduction requires an approved UK Biobank project and locally prepared inputs that satisfy the contracts below.

## Repository map

| File | Manuscript analysis | Main output |
|---|---|---|
| `00_run_matlab_pipeline.m` | Primary MATLAB stages 2-6 | Runs the primary workflow |
| `01_extract_ukb_tabular_data.R` | UKB field extraction | Study-specific CSV files and field manifest |
| `02_compute_sfc_dsfc_glasser.m` | Static and dynamic SFC | Participant x parcel SFC/DSFC MAT file |
| `03_prepare_sfc_analysis_tables.m` | Glasser and Yeo-7 aggregation | MAT and CSV analysis tables |
| `04_run_demographic_age_associations.m` | Age and demographic LMMs | Demographic/age result MAT file |
| `05_run_phenotype_lifestyle_associations.m` | Cognitive, health, and lifestyle LMMs | Association result MAT file |
| `06_run_genetic_associations.m` | PRS, APOE, and AD-variant LMMs | Genetic result MAT file |
| `07_prepare_gwas_inputs.m` | GWAS phenotype preparation | PLINK phenotype TSV |
| `07_run_gwas.sh` | Chromosome-wise GWAS | PLINK 2 result files |
| `08_run_mediation_analysis.R` | Statistical indirect effects | Parameter-estimate CSV |
| `09_run_robustness_analyses.m` | Window-length and SC-measure sensitivity | Robustness summary CSV |

The original exploratory scripts remain recoverable from Git history but are removed from the current tree because they duplicated the cleaned workflow and contained obsolete cluster paths. `utilities/` contains the two mixed-model engines used by the main scripts.

## Software

- MATLAB with Statistics and Machine Learning Toolbox and Parallel Computing Toolbox.
- R with `data.table`, `lmerTest`, and `lavaan`.
- PLINK 2 for GWAS.
- A SLURM cluster is optional; `submit_sfc_dsfc.slurm` is an example, not a requirement.

Run `Rscript R_packages.R` to report whether the required R packages are installed and record `sessionInfo()` when generating the final release. The original files refer to both MATLAB R2018b and MATLAB 2025a; this version discrepancy is listed in `docs/REPRODUCIBILITY_NOTES.md` and should be resolved before publication.

## Expected data layout

The paths can be changed in `config.m`, but the default contract is:

```text
authorized_ukb_data/
├── atlas/
│   └── Yeo7MMP.mat
├── imaging/
│   ├── ukb_glasser_ts_2_0.mat
│   └── ukb_glasser_s1_tract_2_0_streamline_count.mat
├── tabular/
│   ├── UKB_Basic.csv
│   ├── UKB_BrainMRICov.csv
│   ├── UKB_Outcomes_2_0.csv
│   └── UKB_Lifestyle_0_0.csv
└── genetics/
    ├── UKB_GeneCov.csv
    ├── UKB_GeneCov_MRI.txt
    ├── UKB_PRS.tsv
    ├── APOE4_genetype.txt
    └── AD_genetype.txt
```

All participant tables must contain `eid`, or one unambiguous identifier that the relevant script can rename to `eid`.

### Imaging MAT variables

`ukb_glasser_ts_2_0.mat`:

- `fMRI_Glasser_TS_2_0`: N x 1 cell array. Each cell is a 360 x time or time x 360 BOLD matrix.
- `fMRI_Glasser_eID_2_0`: N x 1 participant IDs.

`ukb_glasser_s1_tract_2_0_streamline_count.mat`:

- `connectome_streamline_count_10M_2_0`: M x 1 cell array of structural-connectivity matrices.
- `dMRI_Glasser_S1_eID_2_0`: M x 1 participant IDs.

`Yeo7MMP.mat`:

- `Yeo7MMP`: 360-element vector with integer labels 1-7.

The computation script uses native MATLAB correlations and no longer depends on undocumented local functions such as `DynFunConn`, `SDFC`, `icatb_mat2vec`, or `parcel_to_surface`.

## Quick start

Clone the repository and create a private local configuration:

```bash
cp config_example.m config.m
```

Edit `DATA_ROOT`, optional robustness paths, worker count, and GWAS paths in `config.m`. The file is ignored by Git.

If starting from a legacy UKB bulk table, extract the documented fields:

```bash
export UKB_BULK_FILE=/approved/path/ukb_bulk.csv
export UKB_EXTRACT_DIR=/approved/path/tabular
Rscript 01_extract_ukb_tabular_data.R
```

This extraction preserves the source column names. Harmonize the protected local tables to the analysis-ready names in `docs/DATA_DICTIONARY.md` before running MATLAB; do not upload participant rows.

Run the primary MATLAB workflow from the repository root:

```matlab
run('00_run_matlab_pipeline.m')
```

The default primary analysis uses a 30-TR window and a 5-TR step. To generate 40-TR and 50-TR sensitivity inputs, change `WINDOW_LENGTH` in the private `config.m`, run `02_compute_sfc_dsfc_glasser.m` for each setting, then run `03_prepare_sfc_analysis_tables.m` and `09_run_robustness_analyses.m`.

Run mediation after the 30-TR parcel table exists:

```bash
export UKB_DATA_ROOT=/approved/path/authorized_ukb_data
export UKB_DERIVED_DIR=/path/to/repository/derived
Rscript 08_run_mediation_analysis.R
```

For a fast smoke test of the R script, reduce bootstrap draws without changing the committed default:

```bash
UKB_BOOTSTRAP_SAMPLES=20 Rscript 08_run_mediation_analysis.R
```

## GWAS

First run `07_prepare_gwas_inputs.m`. Then set the paths required by the shell script:

```bash
export PLINK2_BIN=/path/to/plink2
export BFILE_PATTERN='/approved/path/qc/ukb_imp_chr%d'
export PHENO_FILE=/path/to/derived/gwas/UKB_sfc_dsfc_phenotypes_30TR.tsv
export COVAR_FILE=/approved/path/genetics/UKB_GeneCov_MRI.txt
export OUTPUT_DIR=/path/to/derived/gwas/results
bash 07_run_gwas.sh
```

The genotype QC thresholds reported in the manuscript are call rate >=95%, minor-allele frequency >=0.1%, and Hardy-Weinberg equilibrium P >=1e-5. Genotype QC itself is upstream of this repository and must be documented for the authorized UKB release used by the analyst.

## Outputs

Generated data are written under `derived/` by default and are ignored by Git:

```text
derived/
├── ukb_glasser_sfc_dsfc_streamline_count_30TR.mat
├── analysis_tables/
├── association_results/
├── gwas/
├── mediation_results/
└── robustness_results/
```

Each output uses participant IDs only within the authorized environment. Do not commit generated participant-level files.

## Lightweight checks

`python3 tests/check_repository.py` verifies the expected script set, English-only code text, absence of private absolute paths, shell syntax, and balanced MATLAB code blocks. GitHub Actions runs these checks and parses all maintained R scripts without requiring UK Biobank data.

## Analysis conventions

- Participants are matched across imaging modalities using `eid`.
- Functional connectivity is Pearson correlation within each sliding window.
- Regional SFC is Spearman correlation between structural and functional connection profiles after excluding the regional self-connection.
- Zero structural edges and negative functional edges are retained.
- Functional correlations are not Fisher transformed before SFC estimation.
- Participants with mean framewise displacement above 0.20 mm are removed in association scripts.
- Continuous predictors and outcomes are standardized inside the mixed-model utility; binary 0/1 predictors retain their scale.
- Assessment centre is modeled as a random intercept.
- The scripts return raw model P values. Apply the exact family-wise correction described for each manuscript analysis when producing tables and figures.

## Reproducibility notes

Known manuscript/code decisions that require author confirmation are tracked in `docs/REPRODUCIBILITY_NOTES.md`. The most important is whether streamline counts were analyzed without weighting, as stated in the current Methods, or normalized by parcel-pair surface size, as done in the original calculation script. The cleaned default follows the current Methods and records the chosen setting in every SFC/DSFC output MAT file.

## Data availability

UK Biobank data are available through the UK Biobank Access Management System to approved researchers. This repository does not redistribute participant-level imaging, phenotype, lifestyle, or genetic data. The study used UK Biobank application 19542.
