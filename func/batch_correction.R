process_batch_correction <- function() {
  eData <- read_tsv("./output/merged_all_gene_matrix.txt")
  pData <- read_tsv("./output/merged_all_tissue_type_matrix.txt")
  eMatrix <- as.matrix(column_to_rownames(eData, "Symbol"))
  batchList <- pData$batch
  mod <- model.matrix(~as.factor(type), data=pData)
  rownames(mod) <- pData$sample
  bcData <- ComBat(dat=eMatrix, batch=batchList, mod=mod)
  bcMatrix <- cbind(eData[1], as.data.frame(bcData))
  write_tsv(bcMatrix, "./output/all_gene_matrix_after_batch_correction.txt")
  return("./output/all_gene_matrix_after_batch_correction.txt")
}