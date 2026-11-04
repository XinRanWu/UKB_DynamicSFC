# Smoking -> SFC -> fluid intelligence and income -> DSFC -> fluid intelligence
# Primary and socioeconomic-adjusted models use explicit lavaan covariates.

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
source(file.path(root,"R","functions","io.R")); source(file.path(root,"R","functions","lavaan_models.R"))
cfg <- load_project_config(root)

basic <- read_analysis_table(cfg$basic_file,"basic/MRI covariates")
outcomes <- read_analysis_table(cfg$outcome_file,"outcomes")
lifestyle <- read_analysis_table(cfg$lifestyle_file,"lifestyle")
ses <- read_analysis_table(cfg$ses_file,"socioeconomic covariates")
coupling <- read_analysis_table(cfg$coupling_360_file,"parcel coupling")

sfc_parcels <- c(46,67,149,151,247,253,260,329)
dsfc_parcels <- c(33,77,108,283)
sfc_columns <- sprintf("SFC_P%03d",sfc_parcels); dsfc_columns <- sprintf("DSFC_P%03d",dsfc_parcels)
missing <- setdiff(c(sfc_columns,dsfc_columns),names(coupling)); if(length(missing))stop("Missing mediator parcels: ",paste(missing,collapse=", "))
mediators <- data.frame(eid=coupling$eid,
                        M_sfc=rowMeans(coupling[,sfc_columns,drop=FALSE]),
                        M_dsfc=rowMeans(coupling[,dsfc_columns,drop=FALSE]))
d <- inner_join_eid(basic,outcomes,lifestyle,ses,mediators)

smoking_source <- if("SmokingStatus_0_0" %in% names(d)) "SmokingStatus_0_0" else "x20116"
income_source <- if("Income_0_0" %in% names(d)) "Income_0_0" else "x738"
outcome <- if("Fluid_intelligence_score_2_0" %in% names(d)) "Fluid_intelligence_score_2_0" else "x20016"
d$Smoking_ever <- ifelse(is.na(d[[smoking_source]]) | d[[smoking_source]]<0,NA,
                         ifelse(d[[smoking_source]]==0,0,ifelse(d[[smoking_source]]%in%c(1,2),1,NA)))
d$Income_analysis <- ifelse(is.na(d[[income_source]]) | d[[income_source]]<0,NA,as.numeric(d[[income_source]]))

base_covariates <- c("AgeAttend_2_0","Sex_0_0","Race_2","Race_3","Race_4","BMI_2_0","Centre_0_0")
ses_covariates <- c("Education4_0","Occupation_SOC9_0","Employment_retired_0")
brain_covariates <- c("HeadMotion_2_0","Vol_WB_TIV_2_0","SNR_2_0")
categorical <- c("Centre_0_0","Education4_0","Occupation_SOC9_0")
output <- file.path(cfg$derived_root,"results","08_lifestyle_indirect_effects")

# Require SES variables in the shared complete-case samples so attenuation is
# not confounded by a changing sample between paired sensitivity models.
common_smoking <- unique(c("Smoking_ever","M_sfc",outcome,base_covariates,ses_covariates,brain_covariates))
common_income <- unique(c("Income_analysis","M_dsfc",outcome,base_covariates,ses_covariates,brain_covariates))
smoking_data <- d[stats::complete.cases(d[,common_smoking]),]
income_data <- d[stats::complete.cases(d[,common_income]),]

fit_single_mediator(smoking_data,"Smoking_ever","M_sfc",outcome,base_covariates,
                    brain_covariates,categorical,standardize_x=FALSE,bootstrap=cfg$bootstrap_indirect,
                    seed=cfg$seed,analysis_name="smoking_sfc_primary",output_dir=output)
fit_single_mediator(smoking_data,"Smoking_ever","M_sfc",outcome,c(base_covariates,ses_covariates),
                    brain_covariates,categorical,standardize_x=FALSE,bootstrap=cfg$bootstrap_indirect,
                    seed=cfg$seed,analysis_name="smoking_sfc_ses_adjusted",output_dir=output)
fit_single_mediator(income_data,"Income_analysis","M_dsfc",outcome,base_covariates,
                    brain_covariates,categorical,standardize_x=TRUE,bootstrap=cfg$bootstrap_indirect,
                    seed=cfg$seed,analysis_name="income_dsfc_primary",output_dir=output)
fit_single_mediator(income_data,"Income_analysis","M_dsfc",outcome,c(base_covariates,ses_covariates),
                    brain_covariates,categorical,standardize_x=TRUE,bootstrap=cfg$bootstrap_indirect,
                    seed=cfg$seed,analysis_name="income_dsfc_ses_adjusted",output_dir=output)
