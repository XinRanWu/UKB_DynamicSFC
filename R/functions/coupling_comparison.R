fit_age_coupling_comparison <- function(coupling, basic,
                                        networks=c("VN","SMN","DAN","VAN","LN","FPN","DMN")) {
  require_packages(c("dplyr","tidyr","lmerTest","emmeans"))
  sfc <- paste0("SFC_",networks); dsfc <- paste0("DSFC_",networks)
  covariates <- c("Sex_0_0","HeadMotion_2_0","Race_2","Race_3","Race_4",
                  "BMI_2_0","Vol_WB_TIV_2_0","SNR_2_0","Centre_0_0")
  needed <- c("eid","AgeAttend_2_0",covariates)
  d <- inner_join_eid(coupling[,c("eid",sfc,dsfc)],basic[,needed])
  d <- d[stats::complete.cases(d[,c(needed,sfc,dsfc)]),]
  d$Participant <- factor(d$eid); d$Age_z <- zscore(d$AgeAttend_2_0,"age")
  d$Sex_0_0 <- factor(d$Sex_0_0); d$Centre_0_0 <- factor(d$Centre_0_0)
  numeric_covariates <- setdiff(covariates,c("Sex_0_0","Centre_0_0"))
  d[numeric_covariates] <- lapply(d[numeric_covariates],zscore)
  key <- data.frame(measure=c(sfc,dsfc),CouplingType=rep(c("SFC","DSFC"),each=7),
                    Network=rep(networks,2),stringsAsFactors=FALSE)
  long <- tidyr::pivot_longer(d,cols=dplyr::all_of(c(sfc,dsfc)),names_to="measure",values_to="Coupling")
  long <- dplyr::left_join(long,key,by="measure")
  long <- dplyr::group_by(long,.data$CouplingType,.data$Network)
  long <- dplyr::mutate(long,Coupling_z=as.numeric(scale(.data$Coupling)))
  long <- dplyr::ungroup(long)
  long$CouplingType <- factor(long$CouplingType,levels=c("SFC","DSFC")); long$Network <- factor(long$Network,levels=networks)
  cov_rhs <- paste(covariates,collapse=" + ")
  formula <- stats::as.formula(paste("Coupling_z ~ Age_z * CouplingType * Network +",cov_rhs,"+ (1|Participant)"))
  fit <- lmerTest::lmer(formula,data=long,REML=FALSE,
                        control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  omnibus <- as.data.frame(lmerTest::anova(fit,type=3,ddf="Satterthwaite"))
  omnibus$term <- rownames(omnibus); rownames(omnibus) <- NULL
  trends <- emmeans::emtrends(fit,~CouplingType|Network,var="Age_z",lmer.df="asymptotic")
  slopes <- as.data.frame(summary(trends,infer=c(TRUE,TRUE),adjust="none"))
  differences <- as.data.frame(summary(emmeans::contrast(trends,"pairwise",by="Network",adjust="none"),
                                       infer=c(TRUE,TRUE),adjust="none"))
  differences$p_FDR <- stats::p.adjust(differences$p.value,"BH")
  list(fit=fit,omnibus=omnibus,slopes=slopes,differences=differences,n=nrow(d),long=long)
}
