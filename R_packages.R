# R package versions used by the cleaned analysis scripts.
# Install packages manually in an isolated project library, then record the
# exact resolved versions with sessionInfo() for the final release.

required_packages <- c(
  "data.table",  # tabular extraction and joins
  "lmerTest",    # mixed-effects residualization
  "lavaan"       # statistical indirect-effect models
)

installed <- utils::installed.packages()
status <- data.frame(
  package = required_packages,
  installed = required_packages %in% rownames(installed),
  version = vapply(required_packages, function(package_name) {
    if (package_name %in% rownames(installed)) {
      as.character(utils::packageVersion(package_name))
    } else {
      NA_character_
    }
  }, character(1))
)
print(status, row.names = FALSE)
