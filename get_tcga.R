process_tcga <- function(tcgaData) {
  basedir <- getwd()
  symbol <- as.data.frame(tcgaData@rowRanges$external_gene_name)
  colnames(symbol) <- "Symbol"
  exprSet <- bind_cols(symbol, as.data.frame(assay(tcgaData)))
  write_tsv(exprSet, "./GDCdata/exprSet_matrix.txt")
  write_lines(paste(basedir,"/GDCdata/exprSet_matrix.txt", sep=""),"../result.txt", append=TRUE)
  save(tcgaData, file="./GDCdata/tcgaData.rda")
}