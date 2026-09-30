# 08_run_mediation_analysis.R
# Statistical indirect-effect analysis described in the manuscript.
# This script replaces the exploratory legacy notebook with one reproducible
# path for the age/APOE -> SFC or DSFC -> fluid-intelligence models.

required_packages <- c("data.table", "lmerTest", "lavaan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Install the versions documented in R_packages.R."
  )
}

data_root <- Sys.getenv("UKB_DATA_ROOT", unset = "data")
derived_dir <- Sys.getenv("UKB_DERIVED_DIR", unset = "derived")
bootstrap_samples <- as.integer(Sys.getenv("UKB_BOOTSTRAP_SAMPLES", unset = "5000"))

input_files <- c(
  basic = file.path(data_root, "tabular", "UKB_Basic.csv"),
  mri = file.path(data_root, "tabular", "UKB_BrainMRICov.csv"),
  outcomes = file.path(data_root, "tabular", "UKB_Outcomes_2_0.csv"),
  apoe = file.path(data_root, "genetics", "APOE4_genetype.txt"),
  coupling = file.path(derived_dir, "analysis_tables", "brain_data_360p_30TR.csv")
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs: ", paste(missing_files, collapse = ", "))
}

tables <- lapply(input_files, data.table::fread, na.strings = c("", "NA", "NaN"))
for (table_name in names(tables)) {
  if (!("eid" %in% names(tables[[table_name]]))) {
    candidate <- intersect(c("IID", "iid", "ID", "participant_id"), names(tables[[table_name]]))
    if (length(candidate) != 1) {
      stop("Could not identify one participant-ID column in ", input_files[[table_name]])
    }
    data.table::setnames(tables[[table_name]], candidate, "eid")
  }
}

analysis_data <- Reduce(
  function(left, right) merge(left, right, by = "eid", all = FALSE),
  tables
)

fluid_candidates <- c("Fluid_intelligence_score_2_0", "x20016")
fluid_name <- intersect(fluid_candidates, names(analysis_data))[1]
if (is.na(fluid_name)) {
  stop("Fluid-intelligence variable not found; expected x20016 or Fluid_intelligence_score_2_0.")
}
data.table::setnames(analysis_data, fluid_name, "fluid_intelligence")

snr_candidates <- c("SNR_2_0", "SNR_clean_2_0")
snr_name <- intersect(snr_candidates, names(analysis_data))[1]
if (is.na(snr_name)) {
  stop("SNR covariate not found; expected SNR_2_0 or SNR_clean_2_0.")
}

sfc_parcels <- c(4, 5, 6, 13, 142, 152, 183, 193, 332)
dsfc_parcels <- c(38, 45, 74, 81, 100, 115, 116, 179, 335)
mediator_columns <- list(
  SFC = sprintf("SFC_P%03d", sfc_parcels),
  DSFC = sprintf("DSFC_P%03d", dsfc_parcels)
)

required_base <- c(
  "AgeAttend_2_0", "APOE4_dosage", "fluid_intelligence", "Sex_0_0",
  "HeadMotion_2_0", "Race_1", "Race_2", "Race_3", "Race_4", "BMI_2_0",
  "Vol_WB_TIV_2_0", snr_name, "Centre_0_0"
)
missing_columns <- setdiff(
  c(required_base, unlist(mediator_columns, use.names = FALSE)),
  names(analysis_data)
)
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

result_dir <- file.path(derived_dir, "mediation_results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
all_results <- list()

for (mediator_name in names(mediator_columns)) {
  parcel_names <- mediator_columns[[mediator_name]]
  model_data <- analysis_data[, ..required_base]
  model_data[, mediator_raw := rowMeans(analysis_data[, ..parcel_names], na.rm = TRUE)]
  model_data[!is.finite(mediator_raw), mediator_raw := NA_real_]
  model_data <- stats::na.omit(model_data)
  if (nrow(model_data) == 0) {
    stop("No complete cases remain for the ", mediator_name, " mediation model.")
  }

  residual_formula <- stats::as.formula(paste(
    "mediator_raw ~ Sex_0_0 + HeadMotion_2_0 + Race_1 + Race_2 +",
    "Race_3 + Race_4 + BMI_2_0 +", snr_name,
    "+ Vol_WB_TIV_2_0 + (1 | Centre_0_0)"
  ))
  residual_model <- lmerTest::lmer(residual_formula, data = model_data, REML = TRUE)
  model_data[, mediator := as.numeric(scale(stats::residuals(residual_model)))]
  model_data[, AgeAttend_2_0 := as.numeric(scale(AgeAttend_2_0))]
  model_data[, fluid_intelligence := as.numeric(scale(fluid_intelligence))]

  mediation_model <- "
    fluid_intelligence ~ direct_age*AgeAttend_2_0 + direct_apoe*APOE4_dosage + b*mediator
    mediator ~ a_age*AgeAttend_2_0 + a_apoe*APOE4_dosage
    indirect_age := a_age*b
    indirect_apoe := a_apoe*b
    total_age := direct_age + indirect_age
    total_apoe := direct_apoe + indirect_apoe
  "

  fit <- lavaan::sem(
    mediation_model,
    data = model_data,
    se = "bootstrap",
    bootstrap = bootstrap_samples,
    fixed.x = FALSE
  )
  estimates <- lavaan::parameterEstimates(
    fit, standardized = TRUE, ci = TRUE, boot.ci.type = "perc"
  )
  estimates$mediator <- mediator_name
  estimates$n <- nrow(model_data)
  all_results[[mediator_name]] <- estimates
}

combined_results <- data.table::rbindlist(all_results, fill = TRUE)
output_file <- file.path(result_dir, "age_apoe_fluid_intelligence_mediation.csv")
data.table::fwrite(combined_results, output_file)
message("Mediation results saved to: ", normalizePath(output_file))
