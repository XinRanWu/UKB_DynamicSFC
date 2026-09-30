# Analysis-ready data contract

`01_extract_ukb_tabular_data.R` extracts source columns and writes a field manifest. It intentionally preserves the column names of the authorized UKB export. Before running MATLAB association scripts, harmonize those protected local tables to the names below. This step remains local because historical UKB converters and current RAP exports use different column conventions.

## Shared identifier

Every table must contain one row per participant and an `eid` column with the same numeric or string type across files.

## Demographic table

Default file: `tabular/UKB_Basic.csv`

| Analysis column | UKB source | Notes |
|---|---|---|
| `eid` | participant identifier | Join key |
| `AgeAttend_2_0` | Field 21003, imaging instance | Age at imaging |
| `Sex_0_0` | Field 31 | Binary coding must be documented |
| `BMI_2_0` | Field 21001, imaging instance | Continuous |
| `Centre_0_0` | Field 54 at the imaging visit | Random-intercept site; confirm the local instance suffix |
| `Race_1` to `Race_4` | Field 21000 | Study-specific indicator variables; document reference group |

## MRI covariate table

Default file: `tabular/UKB_BrainMRICov.csv`

| Analysis column | UKB source | Notes |
|---|---|---|
| `eid` | participant identifier | Join key |
| `HeadMotion_2_0` | Field 25741 | Exclude values above 0.20 mm |
| `SNR_2_0` or `SNR_clean_2_0` | Field 25744 or locally derived equivalent | Document whether inverse temporal SNR was transformed |
| `Vol_WB_TIV_2_0` | Field 26521 | Total intracranial volume |
| `t1_vol_scaling_2_0` | local UKB imaging covariate | Record the upstream derivation and units |

## Outcome and lifestyle tables

Default files: `tabular/UKB_Outcomes_2_0.csv` and `tabular/UKB_Lifestyle_0_0.csv`.

- Keep `eid` plus analysis variables named `x<FieldID>` unless a main script documents a derived suffix.
- The complete study manifests are `UKB_Outcomes_List.csv` and `UKB_Lifestyle_List.csv`.
- `x1160_L` is a derived long-sleep indicator defined as Field 1160 >=9 hours; the archived code also defined a short-sleep indicator below 6 hours.
- Record recoding, missing-value treatment, direction, and units for every derived binary or ordinal variable in the final release.

## Genetic tables

| File | Required columns |
|---|---|
| `UKB_GeneCov.csv` | `eid`, `x22000_0_0`, and ancestry PCs `x22009_0_1` to `x22009_0_10` |
| `UKB_PRS.tsv` | `eid` plus the tested PRS columns |
| `APOE4_genetype.txt` | `eid`, `APOE4_dosage` |
| `AD_genetype.txt` | `eid` plus candidate AD-variant dosage columns |

The final public release should include a de-identified schema-only example for each table. Do not include participant rows.
