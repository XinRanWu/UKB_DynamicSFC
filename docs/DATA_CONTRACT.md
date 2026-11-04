# Analysis ready data contract

Participant-level UK Biobank data cannot be distributed. Reproduction therefore starts from authorized local files with the schemas below. Every table must contain one unique, nonmissing `eid`; scripts stop on duplicate identifiers rather than silently multiplying rows.

## Imaging MAT files

- Glasser fMRI: participant IDs plus an N by 1 cell array. Each cell contains a 360 by time matrix or its transpose.
- Glasser SC: participant IDs plus an N by 1 cell array of 360 by 360 symmetric matrices.
- GSR fMRI and mean-FA SC use the same participant-ID and cell-array contracts.
- `Yeo7MMP.mat` contains one integer network label from 1 to 7 for each parcel.

Exact MAT variable names are set in `config/config.m`; they are never inferred.

## Basic and MRI covariates

Required columns are `eid`, `AgeAttend_2_0`, `Sex_0_0`, `Race_1` to `Race_4`, `BMI_2_0`, `HeadMotion_2_0`, `SNR_2_0`, `Vol_WB_TIV_2_0`, and `Centre_0_0`. The code excludes mean framewise displacement above 0.20 mm.

## Phenotype and lifestyle tables

The outcome and lifestyle tables contain `eid` followed by the variables to be tested. The exact final model columns are fixed in `metadata/phenotypes_used.csv` (75 phenotypes) and `metadata/lifestyle_used.csv` (56 environmental/lifestyle measures); the scripts fail if any are absent. Descriptions and categories are provided in `UKB_Outcomes_List.csv`, `UKB_Lifestyle_List.csv`, and Supplementary Table S1. Negative UKB response codes must be converted to missing values during local harmonization.

## Socioeconomic covariates

`UKB_SES_encoded_0.csv` contains `Education4_0`, `Occupation_SOC9_0`, and `Employment_retired_0`. The coding used in the revised manuscript is documented in the Supplementary Methods.

## Genetic inputs

- PRS: `eid` plus the 90 tested PRS variables; AD-PRS is identified by `26206` or `alzheimer` in its column name.
- Genetic covariates: `eid` and `x22009_0_1` to `x22009_0_10`.
- APOE: `eid` and `APOE4_dosage` coded 0, 1, or 2.

The genotyping array variable `x22000_0_0` is not used in the final PRS/SEM scripts because the final reported analysis excluded this anomalous variable. Genetic ancestry is controlled using the first ten principal components.

## Local harmonization

`R/00_extract_ukb_fields.R` can extract field columns from a legacy UKB bulk export, but it deliberately does not guess local converter-specific names or response-code transformations. Harmonize extracted columns to the exact names above and record that mapping inside the authorized analysis environment.
