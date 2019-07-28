process_ciber <- function(geneMatrixPath) {
  source("./func/process_ciber.R")
  do_cibersort(geneMatrixPath)
}

do_cibersort <- function (filePath){
  source("./func/supportFunc_cibersort.R")
  results <- CIBERSORT("./data/lm22.txt",filePath,perm=1000, QN=TRUE)
  save(results, file="./output/cibersort.rda")
  tResults <- as.data.frame(t(results))
  result_to_save <- rownames_to_column(tResults, "immCell")
  write_tsv(result_to_save,"./output/cibersort_all_result.txt")
}