process_mcp <- function(geenMatrixPath) {
  eData <- read_tsv(geenMatrixPath)
  eMatrix <- as.matrix(column_to_rownames(eData, "Symbol"))
  res <- MCPcounter::MCPcounter.estimate(eMatrix, featuresType = c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[2])
  res_to_save <- as.data.frame(rownames_to_column(as.data.frame(res), "Type"))
  write_tsv(res_to_save, "./output/mcp_all_result.txt")
}