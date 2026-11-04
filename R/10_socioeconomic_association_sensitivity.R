# Network-level smoking and income associations after additional adjustment
# for education, occupational class, and retirement status.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
source(file.path(root,"R","functions","io.R")); cfg <- load_project_config(root)
require_packages(c("lmerTest"))
basic <- read_analysis_table(cfg$basic_file,"basic/MRI covariates")
lifestyle <- read_analysis_table(cfg$lifestyle_file,"lifestyle")
ses <- read_analysis_table(cfg$ses_file,"socioeconomic covariates")
coupling <- read_analysis_table(cfg$coupling_7n_file,"network coupling")
d <- inner_join_eid(basic,lifestyle,ses,coupling)
smoking <- if("SmokingStatus_0_0"%in%names(d))"SmokingStatus_0_0"else"x20116"
income <- if("Income_0_0"%in%names(d))"Income_0_0"else"x738"
d$Smoking_ever <- ifelse(is.na(d[[smoking]])|d[[smoking]]<0,NA,ifelse(d[[smoking]]==0,0,1))
d$Income_analysis <- ifelse(is.na(d[[income]])|d[[income]]<0,NA,as.numeric(d[[income]]))
d$AgeSq <- zscore(d$AgeAttend_2_0)^2
d$Sex_0_0 <- factor(d$Sex_0_0); d$Centre_0_0 <- factor(d$Centre_0_0)
d$Education4_0 <- factor(d$Education4_0); d$Occupation_SOC9_0 <- factor(d$Occupation_SOC9_0)
outcomes <- c(paste0("SFC_",c("VN","SMN","DAN","VAN","LN","FPN","DMN")),
              paste0("DSFC_",c("VN","SMN","DAN","VAN","LN","FPN","DMN")))
covariates <- c("AgeAttend_2_0","AgeSq","Sex_0_0","Race_2","Race_3","Race_4","BMI_2_0",
                "HeadMotion_2_0","SNR_2_0","Vol_WB_TIV_2_0","Education4_0",
                "Occupation_SOC9_0","Employment_retired_0")
rows <- list(); index <- 0L
for(x in c("Smoking_ever","Income_analysis")) for(y in outcomes) {
  variables <- c(x,y,covariates,"Centre_0_0"); model_data <- d[stats::complete.cases(d[,variables]),variables]
  model_data[[y]] <- zscore(model_data[[y]],y)
  if(x=="Income_analysis") model_data[[x]] <- zscore(model_data[[x]],x)
  formula <- stats::as.formula(paste(y,"~",paste(c(x,covariates),collapse=" + "),"+ (1|Centre_0_0)"))
  fit <- lmerTest::lmer(formula,data=model_data,REML=TRUE)
  coefficient <- summary(fit)$coefficients[x,]
  ci <- stats::confint(fit,parm=x,method="Wald")
  index <- index+1L
  rows[[index]] <- data.frame(predictor=x,outcome=y,N=nrow(model_data),beta=coefficient["Estimate"],
                              SE=coefficient["Std. Error"],DF=coefficient["df"],t=coefficient["t value"],
                              p=coefficient["Pr(>|t|)"],CI_lower=ci[1],CI_upper=ci[2])
}
results <- do.call(rbind,rows)
results$coupling_type <- ifelse(grepl("^SFC_",results$outcome),"SFC","DSFC")
results$p_FDR <- ave(results$p,results$predictor,results$coupling_type,FUN=function(p)stats::p.adjust(p,"BH"))
output <- file.path(cfg$derived_root,"results","10_socioeconomic_sensitivity")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(results,file.path(output,"network_associations.csv"),row.names=FALSE)
