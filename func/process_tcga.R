process_tcga <- function() {
  #tcga_project <- "TCGA-STAD"
  #tcga_experimental_strategy <- "RNA-Seq"
  #tcga_data_category <- "Transcriptome Profiling"
  #tcga_data_type <- "Gene Expression Quantification"
  #tcga_workflow_type <- "HTSeq - FPKM"
  
  #tcga_experimental_strategy <- "RNA-Seq"
  #tcga_data_category <- "Transcriptome Profiling"
  #tcga_data_type <- "Gene Expression Quantification"
  #tcga_workflow_type <- "HTSeq - FPKM"
  #tcga_project <- c("TCGA-LUAD", "TCGA-LUSC")
  
  
  setwd("./workspace")
  query <- GDCquery(project = tcga_project, 
                    data.category = tcga_data_category, 
                    data.type = tcga_data_type, 
                    workflow.type = tcga_workflow_type, 
                    experimental.strategy = tcga_experimental_strategy)
  GDCdownload(query)
  tcgaData <- GDCprepare(query)
  complete_tcga(tcgaData)
  setwd("..")
}

fpkmToTpm <- function(x)
{
  lDenom <- log(sum(x))
  exp(log(x) - lDenom + log(1e6))
}

complete_tcga <- function(tcgaData) {
  basedir <- getwd()
  symbol <- as.data.frame(tcgaData@rowRanges$external_gene_name)
  colnames(symbol) <- "Symbol"
  exprSet <- bind_cols(symbol, as.data.frame(assay(tcgaData)))
  # convert fpkm to tpm
  exprSet <- cbind(symbol, apply(exprSet[,-1], 2, fpkmToTpm))
  exprSet_noNA_nodup <- exprSet %>% group_by_at(1) %>% summarise_all(funs(mean))
  
  write_tsv(exprSet, "./GDCdata/exprSet_matrix.txt")
  write_tsv(exprSet_noNA_nodup, "./GDCdata/exprSet_noNA_noDup_matrix.txt")
  exprSet_log2tpm <- log2(column_to_rownames(exprSet_noNA_nodup, "Symbol")+1)
  exprSet_log2tpm <- rownames_to_column(exprSet_log2tpm, "Symbol")
  write_tsv(exprSet_log2tpm, "./GDCdata/exprSet_log2tpm_final.txt")
  write_lines(paste(basedir,"/GDCdata/exprSet_log2tpm_final.txt", sep=""),"../result.txt", append=TRUE)
  save(tcgaData, file="./GDCdata/tcgaData.rda")
}