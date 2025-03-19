library(data.table)
file.vec <- c("db-problems.tsv", "db-loss.tsv", "loss.timing.tsv")
data.list <- list()
for(f in file.vec){
  if(!file.exists(f)){
    options(timeout=600)
    download.file(
      paste0(
        "https://rcdata.nau.edu/genomic-ml/fullpath/",
        f),
      f, cacheOK = TRUE)
  }
  data.list[[f]] <- fread(f)
}
sum(data.list$loss.timing.tsv$seconds, na.rm = TRUE)/60/60/24/365


data.list[["db-loss.tsv"]][
, .(models=.N), by=prob.id
][
  data.list[["db-problems.tsv"]],
  on="prob.id"
]
