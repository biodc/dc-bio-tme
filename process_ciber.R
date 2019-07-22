do_cibersort <- function (filePath){
  source("../func/supportFunc_cibersort.R")
  results <- CIBERSORT(filePath,"../data/lm22.txt", perm=1000, QN=TRUE)
  t <- strsplit(filePath, "/")[[1]]
  fn <- t[length(t)]
  immCell <- rownames(results)
  c <- data.frame(immCell)
  result_to_save <- bind_cols(c, as.data.frame(results))
  write_tsv(result_to_save, path=paste("./", fn, "_ciber_matrix.txt",sep=""))
}