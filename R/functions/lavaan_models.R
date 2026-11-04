encode_covariates <- function(data, covariates, categorical) {
  map <- setNames(vector("list", length(covariates)), covariates)
  for (variable in covariates) {
    if (variable %in% categorical || is.factor(data[[variable]]) || is.character(data[[variable]])) {
      f <- droplevels(factor(data[[variable]]))
      if (nlevels(f) < 2L) { map[[variable]] <- character(); next }
      mm <- stats::model.matrix(~ f)[, -1, drop = FALSE]
      names_new <- make.unique(make.names(paste0("D_", variable, "_", levels(f)[-1])))
      colnames(mm) <- names_new
      data[names_new] <- mm
      map[[variable]] <- names_new
    } else {
      observed <- unique(data[[variable]][is.finite(data[[variable]])])
      if (length(observed) > 2L) data[[variable]] <- zscore(data[[variable]], variable)
      map[[variable]] <- variable
    }
  }
  terms <- unname(unlist(map, use.names = FALSE))
  varying <- terms[vapply(terms, function(x) {
    s <- stats::sd(data[[x]], na.rm = TRUE); is.finite(s) && s > 0
  }, logical(1))]
  list(data = data, map = map, varying = varying)
}

expand_terms <- function(names_in, encoded) {
  intersect(unname(unlist(encoded$map[names_in], use.names = FALSE)), encoded$varying)
}

fit_single_mediator <- function(data, x, mediator, outcome, general_covariates,
                                mediator_covariates = character(), categorical = character(),
                                standardize_x = TRUE, bootstrap = 2000L, seed = 20261104L,
                                analysis_name = "indirect_effect", output_dir = NULL) {
  require_packages("lavaan")
  variables <- unique(c(x, mediator, outcome, general_covariates, mediator_covariates))
  missing <- setdiff(variables, names(data)); if (length(missing)) stop("Missing variables: ", paste(missing, collapse = ", "))
  d <- data[stats::complete.cases(data[, variables, drop = FALSE]), variables, drop = FALSE]
  if (!nrow(d)) stop("No complete cases for ", analysis_name)
  d$X_model <- if (standardize_x) zscore(d[[x]], x) else as.numeric(d[[x]])
  d$M_model <- zscore(d[[mediator]], mediator)
  d$Y_model <- zscore(d[[outcome]], outcome)
  encoded <- encode_covariates(d, unique(c(general_covariates, mediator_covariates)), categorical)
  d <- encoded$data
  general_terms <- expand_terms(general_covariates, encoded)
  mediator_terms <- unique(c(general_terms, expand_terms(mediator_covariates, encoded)))
  syntax <- paste0(
    "M_model ~ a*X_model", if (length(mediator_terms)) paste0(" + ", paste(mediator_terms, collapse = " + ")) else "", "\n",
    "Y_model ~ c_prime*X_model + b*M_model", if (length(general_terms)) paste0(" + ", paste(general_terms, collapse = " + ")) else "", "\n",
    "indirect := a*b\n", "total := c_prime + indirect\n"
  )
  set.seed(seed)
  fit <- lavaan::sem(syntax, data = d, se = "bootstrap", bootstrap = bootstrap,
                     fixed.x = TRUE, meanstructure = TRUE)
  parameters <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE, boot.ci.type = "perc")
  target <- parameters[parameters$label %in% c("a","b","c_prime") |
                         parameters$lhs %in% c("indirect","total"), ]
  target$analysis <- analysis_name; target$n <- nrow(d)
  if (!is.null(output_dir)) export_lavaan_fit(fit, syntax, parameters, target, output_dir, analysis_name)
  invisible(list(fit = fit, syntax = syntax, parameters = parameters, target = target, n = nrow(d), data = d))
}

fit_age_adprs_sem <- function(data, mediator, outcome, age, ad_prs,
                              general_covariates, mediator_covariates = character(),
                              categorical = character(), bootstrap = 1000L,
                              seed = 20261104L, analysis_name = "age_adprs_sem", output_dir = NULL) {
  require_packages("lavaan")
  variables <- unique(c(mediator,outcome,age,ad_prs,general_covariates,mediator_covariates))
  missing <- setdiff(variables,names(data)); if (length(missing)) stop("Missing variables: ",paste(missing,collapse=", "))
  d <- data[stats::complete.cases(data[,variables,drop=FALSE]),variables,drop=FALSE]
  d$M_model <- zscore(d[[mediator]],mediator); d$Y_model <- zscore(d[[outcome]],outcome)
  d$Age_model <- zscore(d[[age]],age); d$ADPRS_model <- zscore(d[[ad_prs]],ad_prs)
  encoded <- encode_covariates(d,unique(c(general_covariates,mediator_covariates)),categorical); d <- encoded$data
  g <- expand_terms(general_covariates,encoded); m <- unique(c(g,expand_terms(mediator_covariates,encoded)))
  syntax <- paste0(
    "M_model ~ a_age*Age_model + a_adprs*ADPRS_model",if(length(m))paste0(" + ",paste(m,collapse=" + "))else"","\n",
    "Y_model ~ c_age*Age_model + c_adprs*ADPRS_model + b*M_model",if(length(g))paste0(" + ",paste(g,collapse=" + "))else"","\n",
    "indirect_age := a_age*b\nindirect_adprs := a_adprs*b\n",
    "total_age := c_age + indirect_age\ntotal_adprs := c_adprs + indirect_adprs\n")
  set.seed(seed)
  fit <- lavaan::sem(syntax,data=d,se="bootstrap",bootstrap=bootstrap,fixed.x=TRUE,meanstructure=TRUE)
  parameters <- lavaan::parameterEstimates(fit,standardized=TRUE,ci=TRUE,boot.ci.type="perc")
  target_labels <- c("a_age","a_adprs","b","c_age","c_adprs")
  target <- parameters[parameters$label %in% target_labels |
                         parameters$lhs %in% c("indirect_age","indirect_adprs","total_age","total_adprs"),]
  target$analysis <- analysis_name; target$n <- nrow(d)
  if(!is.null(output_dir)) export_lavaan_fit(fit,syntax,parameters,target,output_dir,analysis_name)
  invisible(list(fit=fit,syntax=syntax,parameters=parameters,target=target,n=nrow(d),data=d))
}

export_lavaan_fit <- function(fit, syntax, parameters, target, output_dir, analysis_name) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(parameters,file.path(output_dir,paste0(analysis_name,"_parameters.csv")),row.names=FALSE)
  utils::write.csv(target,file.path(output_dir,paste0(analysis_name,"_target_effects.csv")),row.names=FALSE)
  writeLines(syntax,file.path(output_dir,paste0(analysis_name,"_model.txt")))
  writeLines(capture.output(summary(fit,standardized=TRUE,ci=TRUE,rsquare=TRUE)),
             file.path(output_dir,paste0(analysis_name,"_summary.txt")))
}
