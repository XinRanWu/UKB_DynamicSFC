# Covariate-adjusted network age curves used in Figure 2C.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
source(file.path(root,"R","functions","io.R")); cfg <- load_project_config(root)
require_packages(c("mgcv"))
basic <- read_analysis_table(cfg$basic_file,"basic/MRI covariates")
coupling <- read_analysis_table(cfg$coupling_7n_file,"network coupling")
d <- inner_join_eid(basic,coupling); d <- d[d$HeadMotion_2_0<=0.20,]
d$Sex_0_0 <- factor(d$Sex_0_0); d$Centre_0_0 <- factor(d$Centre_0_0)
outcomes <- c(paste0("SFC_",c("VN","SMN","DAN","VAN","LN","FPN","DMN")),
              paste0("DSFC_",c("VN","SMN","DAN","VAN","LN","FPN","DMN")))
age_grid <- seq(min(d$AgeAttend_2_0,na.rm=TRUE),max(d$AgeAttend_2_0,na.rm=TRUE),length.out=200)
rows <- list(); fits <- list()
for(y in outcomes){
  variables <- c(y,"AgeAttend_2_0","Sex_0_0","HeadMotion_2_0","Centre_0_0")
  model_data <- d[stats::complete.cases(d[,variables]),variables]
  model <- mgcv::gamm(stats::as.formula(paste(y,"~ s(AgeAttend_2_0,k=10) + Sex_0_0 + HeadMotion_2_0")),
                      random=list(Centre_0_0=~1),data=model_data,method="REML")
  fits[[y]] <- model
  for(sex in levels(model_data$Sex_0_0)){
    new_data <- data.frame(AgeAttend_2_0=age_grid,Sex_0_0=factor(sex,levels=levels(model_data$Sex_0_0)),
                           HeadMotion_2_0=mean(model_data$HeadMotion_2_0))
    prediction <- stats::predict(model$gam,newdata=new_data,se.fit=TRUE)
    rows[[length(rows)+1L]] <- data.frame(outcome=y,sex=sex,age=age_grid,fit=prediction$fit,
      SE=prediction$se.fit,CI_lower=prediction$fit-1.96*prediction$se.fit,
      CI_upper=prediction$fit+1.96*prediction$se.fit)
  }
}
output <- file.path(cfg$derived_root,"results","03_age"); dir.create(output,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(do.call(rbind,rows),file.path(output,"gamm_age_curves.csv"),row.names=FALSE)
saveRDS(fits,file.path(output,"gamm_models.rds"))
