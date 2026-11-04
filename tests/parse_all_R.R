files <- list.files(c("R","tests"),pattern="\\.R$",recursive=TRUE,full.names=TRUE)
files <- setdiff(files,"tests/parse_all_R.R")
for(file in files){message("Parsing ",file); parse(file=file)}
message("All R files parsed successfully.")
