# Record the R environment and repository revision beside final outputs.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))), ".."))
source(file.path(root,"R","functions","io.R")); cfg <- load_project_config(root)
output <- file.path(cfg$derived_root,"run_metadata"); dir.create(output,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(utils::sessionInfo()),file.path(output,"R_sessionInfo.txt"))
revision <- tryCatch(system2("git",c("-C",shQuote(root),"rev-parse","HEAD"),stdout=TRUE,stderr=TRUE),
                     error=function(e)paste("Unavailable:",conditionMessage(e)))
writeLines(revision,file.path(output,"git_commit.txt"))
message("Run metadata written to ",output)
