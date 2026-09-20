# UKB_DynamicSFC
This repository contains the data analysis code used for the paper entitled “Stable brain architecture but flexible dynamics: distinct structure–function coupling signatures link aging and genetic Alzheimer's risk to cognitive decline in 34,067 adults.”

The original analyses were conducted using UK Biobank data downloaded underthe legacy UKB data-access framework, prior to the transition to the UK Biobank Research Analysis Platform (RAP).

The script `01_extract_ukb_tabular_data.R` documents the UK Biobank field IDs used in this study and reproduces the extraction of study-specific variables from an authorized legacy UKB bulk tabular export. Participant-level UK Biobank data are not distributed with this repository.

Users working within the current UKB RAP environment should retrieve the corresponding fields from their approved UKB project and map them to the documented variables before running the downstream analysis scripts.
