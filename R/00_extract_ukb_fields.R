# Extract the documented UKB fields from an authorized legacy bulk export.
# This script does not distribute data and does not replace upstream UKB/RAP
# imaging or genotype processing.
requireNamespace("data.table",quietly=TRUE) || stop("Install data.table.")
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
bulk_file <- Sys.getenv("UKB_BULK_FILE",unset=NA_character_)
output_dir <- Sys.getenv("UKB_EXTRACT_DIR",unset=file.path(root,"derived","extracted_fields"))
if(is.na(bulk_file)||!file.exists(bulk_file))stop("Set UKB_BULK_FILE to an authorized UKB export.")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
outcome_manifest <- data.table::fread(file.path(root,"UKB_Outcomes_List.csv"))
lifestyle_manifest <- data.table::fread(file.path(root,"UKB_Lifestyle_List.csv"))
groups <- list(
  basic=c(31,54,21003,21001,21000,26521,25741,25744),
  outcomes=unique(outcome_manifest[["Field ID"]]),
  lifestyle=unique(lifestyle_manifest[["Field ID"]]),
  socioeconomic=c(6138,20277,6142,20119),
  genetic_covariates=c(22009)
)
matches_field <- function(column,id) any(vapply(c(paste0("^",id,"[-.]"),paste0("^f\\.",id,"\\."),
  paste0("^x",id,"_"),paste0("^",id,"$")),grepl,logical(1),x=column))
ukb <- data.table::fread(bulk_file,showProgress=TRUE)
if(!"eid"%in%names(ukb))stop("Rename the participant identifier to eid before extraction.")
for(group in names(groups)){
  keep <- c("eid",names(ukb)[vapply(names(ukb),function(column)
    any(vapply(groups[[group]],function(id)matches_field(column,id),logical(1))),logical(1))])
  data.table::fwrite(ukb[,..keep],file.path(output_dir,paste0("UKB_",group,"_fields.csv")),na="NA")
}
manifest <- data.table::rbindlist(lapply(names(groups),function(group)
  data.table::data.table(group=group,field_id=groups[[group]])))
data.table::fwrite(manifest,file.path(output_dir,"field_manifest.csv"))
