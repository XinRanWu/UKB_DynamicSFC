repo_root <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", grep("^--file=", command, value = TRUE))
  if (length(script) == 1L) return(normalizePath(file.path(dirname(script), "..")))
  normalizePath(".")
}

load_project_config <- function(root = repo_root()) {
  path <- file.path(root, "config", "config.R")
  if (!file.exists(path)) stop("Missing config/config.R; copy config/config.example.R first.")
  environment <- new.env(parent = baseenv())
  sys.source(path, envir = environment)
  if (!exists("config", envir = environment, inherits = FALSE)) stop("config.R must define config.")
  environment$config
}

read_analysis_table <- function(path, name = basename(path), sep = NULL) {
  if (!file.exists(path)) stop("Missing ", name, ": ", path)
  if (is.null(sep)) sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  out <- utils::read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                           stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN"))
  if (!"eid" %in% names(out)) stop(name, " must contain eid.")
  if (anyNA(out$eid) || anyDuplicated(out$eid)) stop(name, " must contain unique nonmissing eid values.")
  out$eid <- as.character(out$eid)
  out
}

inner_join_eid <- function(...) {
  tables <- list(...)
  Reduce(function(x, y) {
    duplicates <- setdiff(intersect(names(x), names(y)), "eid")
    if (length(duplicates)) stop("Duplicate variables before merge: ", paste(duplicates, collapse = ", "))
    merge(x, y, by = "eid", all = FALSE, sort = FALSE)
  }, tables)
}

zscore <- function(x, name = deparse(substitute(x))) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) stop(name, " has zero or invalid variance.")
  as.numeric(scale(x))
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
}
