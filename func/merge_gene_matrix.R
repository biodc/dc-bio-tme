mergeGeneMatrix <- function(targetFilePath) {
  targetFileList <- read_lines(targetFilePath)
  first <- TRUE
  batch <- 1
  for(targetFile in targetFileList) {
    if(str_detect(targetFile, "GDCdata")) {
      tcga_gene_matrix <- read_tsv(targetFile)
      sample <- colnames(tcga_gene_matrix)[-1]
      type <- data.frame(type=sapply(sample, detectTCGA), row.names=NULL)
      tcga_pdata <- data.frame(sample, type, batch)
      batch <- batch+1
    } else {
      pdataPath <- str_replace(targetFile, "matrix_expr_symbol_noNAnoDup", "pdata")
      pdata <- read_tsv(pdataPath)
      sample <- pdata$geo_accession
      tissue <- pdata$source_name_ch1
      type <- data.frame(type=sapply(tissue, detectTumor))
      if(first) {
        finalMatrix <- read_tsv(targetFile)
        colnames(finalMatrix) <- c(colnames(finalMatrix)[1], pdata$geo_accession)
        finalPdata <- data.frame(sample, type, batch)
        batch <- batch+1
        first <- FALSE
      } else {
        temp_matrix <- read_tsv(targetFile)
        colnames(temp_matrix) <- c(colnames(temp_matrix)[1], pdata$geo_accession)
        finalMatrix %<>% inner_join(temp_matrix, by="Symbol")
        finalPdata %<>% union_all(data.frame(sample, type, batch))
        batch <- batch+1
      }
    }
  }
  if(exists("tcga_gene_matrix")){
    finalMatrix %<>% inner_join(tcga_gene_matrix, by="Symbol")
    finalPdata %<>% union_all(tcga_pdata)
  }
  write_tsv(finalMatrix,"./output/merged_all_gene_matrix.txt")
  write_tsv(finalPdata, "./output/merged_all_tissue_type_matrix.txt")
  return(paste(getwd(), "/output/merged_all_gene_matrix.txt", sep=""))
}

detectTumor <- function(x) {
  if(str_detect(x, "normal") || str_detect(x, "adjacent")) {
    return("normal")
  } else {
    return("tumor")
  }
}

detectTCGA <- function(x) {
  targetCode <- str_split(x, "-", simplify = TRUE)[4]
  detectionCode <- str_sub(targetCode, 1, 1)
  if(detectionCode == "0") {
    return("tumor")
  } else {
    return("normal")
  }
}