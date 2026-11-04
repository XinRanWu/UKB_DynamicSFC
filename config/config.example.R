# Copy to config/config.R and edit only local paths and optional run settings.
config <- list(
  data_root = "/path/to/authorized/UKB/data",
  derived_root = file.path(getwd(), "derived"),
  basic_file = "/path/to/authorized/UKB/data/tabular/UKB_BasicMRI.csv",
  outcome_file = "/path/to/authorized/UKB/data/tabular/UKB_Outcomes_2_0.csv",
  lifestyle_file = "/path/to/authorized/UKB/data/tabular/UKB_Lifestyle_0_0.csv",
  ses_file = "/path/to/authorized/UKB/data/tabular/UKB_SES_encoded_0.csv",
  prs_file = "/path/to/authorized/UKB/data/genetics/UKB_PRS.tsv",
  gene_covariate_file = "/path/to/authorized/UKB/data/genetics/UKB_GeneCov.csv",
  coupling_360_file = file.path(getwd(), "derived", "tables", "coupling_360.csv"),
  coupling_7n_file = file.path(getwd(), "derived", "tables", "coupling_7n.csv"),
  bootstrap_indirect = 2000L,
  bootstrap_sem = 1000L,
  seed = 20261104L
)
