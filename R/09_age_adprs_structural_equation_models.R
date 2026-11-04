# Multiple-predictor, single-mediator models:
# age + AD-PRS -> SFC/DSFC -> fluid intelligence.

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
source(file.path(root,"R","functions","io.R")); source(file.path(root,"R","functions","lavaan_models.R"))
cfg <- load_project_config(root)
basic <- read_analysis_table(cfg$basic_file,"basic/MRI covariates")
outcomes <- read_analysis_table(cfg$outcome_file,"outcomes")
prs <- read_analysis_table(cfg$prs_file,"PRS")
gene <- read_analysis_table(cfg$gene_covariate_file,"genetic covariates")
coupling <- read_analysis_table(cfg$coupling_360_file,"parcel coupling")

sfc_columns <- sprintf("SFC_P%03d",c(4,5,6,13,142,152,183,332))
dsfc_columns <- sprintf("DSFC_P%03d",c(38,47,57,59,60,109,179))
mediators <- data.frame(eid=coupling$eid,M_sfc=rowMeans(coupling[,sfc_columns,drop=FALSE]),
                        M_dsfc=rowMeans(coupling[,dsfc_columns,drop=FALSE]))
d <- inner_join_eid(basic,outcomes,prs,gene,mediators)
outcome <- if("Fluid_intelligence_score_2_0" %in% names(d)) "Fluid_intelligence_score_2_0" else "x20016"
ad_candidates <- names(prs)[grepl("26206|alzheimer",names(prs),ignore.case=TRUE)]
if(length(ad_candidates)!=1L)stop("Expected one AD-PRS column; found: ",paste(ad_candidates,collapse=", "))
pcs <- paste0("x22009_0_",1:10)
general <- c("Sex_0_0","BMI_2_0","Centre_0_0",pcs)
brain <- c("HeadMotion_2_0","Vol_WB_TIV_2_0","SNR_2_0")
output <- file.path(cfg$derived_root,"results","09_age_adprs_sem")

# Use one complete-case sample for the paired SFC and DSFC models.
needed <- unique(c("M_sfc","M_dsfc",outcome,"AgeAttend_2_0",ad_candidates,general,brain))
d <- d[stats::complete.cases(d[,needed]),]
fit_age_adprs_sem(d,"M_sfc",outcome,"AgeAttend_2_0",ad_candidates,general,brain,"Centre_0_0",
                  cfg$bootstrap_sem,cfg$seed,"age_adprs_sfc",output)
fit_age_adprs_sem(d,"M_dsfc",outcome,"AgeAttend_2_0",ad_candidates,general,brain,"Centre_0_0",
                  cfg$bootstrap_sem,cfg$seed,"age_adprs_dsfc",output)
