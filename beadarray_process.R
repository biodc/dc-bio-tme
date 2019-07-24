process_geo_ill <- function(gse) {
  basedir <- getwd()
  # if there is non-normalized
  nonNormal <- getGEOSuppFiles(gse, filter_regex = ".*?normalized.*?")
  setwd(gse)
  matrixData <- getGEO(gse, AnnotGPL = FALSE, getGPL = FALSE, destdir = ".")
  fileToName <- gsub(".gz","",attr(matrixData, "name"))
  gplSerial <- matrixData[[1]]@annotation
  gplFiles <- getGEO(gplSerial, destdir = ".")
  gpl <- Table(gplFiles)
  if(length(rownames(nonNormal))>0) {
    gunzip(rownames(nonNormal))
    illData <- read_tsv(gsub(".gz", "", rownames(nonNormal)))
    pvalColumns <- str_detect(colnames(illData), "Detection Pval") %>% which
    sampleColumns <- colnames(illData[,pvalColumns-1])
    sampleTarget <- str_split(colnames(illData[,pvalColumns-1]), pattern = "\\W")[[1]][1]
    probeIDTarget <- colnames(illData[,1])
    rawLumi <- read.ilmn(gsub(".gz", "", rownames(nonNormal)), probeid = probeIDTarget, expr = sampleTarget, other.columns = "Detection")
    exprSet <- as.data.frame(neqc(rawLumi))
  } else {
    rawData <- getGEOSuppFiles(gse, filter_regex = ".*?RAW.*?")
    untar(rownames(rawData))
    sapply(dir(pattern="gz$"), gunzip)
    illData <- read.idat(idatfiles = dir(pattern = "idat$"), bgxfile=dir(pattern="bgx$"))
    illData.norm <- neqc(illData)
    exprSet <- as.data.frame(illData.norm)
    row.names(exprSet) <- (exprSet[1][[1]])
    # mask array and probe
    exprSet <- exprSet[,-c(colnames(exprSet) %>% str_detect("Array") %>% which,colnames(exprSet) %>% str_detect("Probe") %>% which)] 
  }
  probeSet <- rownames(exprSet)
  exprSet$ID <- probeSet
  if(all(c("Gene Symbol", "ID") %in% colnames(gpl))){
    map_probe_symbol <- gpl[,c("ID", "Gene Symbol")]
  } else if (all(c("ORF", "ID") %in% colnames(gpl))){
    map_probe_symbol <- gpl[,c("ID", "ORF")] 
  } else if(all(c("Symbol", "ID") %in% colnames(gpl))) {
    map_probe_symbol <- gpl[,c("ID", "Symbol")]
  } else {
    ann_from_raw <- paste("illumina",illData@annotation, sep = "")
    ann <- annPkgName(ann_from_raw, type="db")
    if(!require(ann, character.only = TRUE)){
      print(paste("Installing ", ann))
      BiocManager::install(ann)
      require(ann, character.only = TRUE)
    }
    map_probe_symbol <- as.data.frame(annotate::getSYMBOL(probeSet, ann_from_raw))
    colnames(map_probe_symbol) <- c("Symbol")
    map_probe_symbol$ID <- rownames(map_probe_symbol)
  }
  if(exists("map_probe_symbol")) {
    exprSet_with_symbol <- merge(map_probe_symbol, exprSet, by="ID")
    final_exprSet <- subset(exprSet_with_symbol, select=-ID)
    colnames(final_exprSet) <- c("Symbol", colnames(final_exprSet[-1]))
    setwd("..")
    if(dir.exists("exprSet_ann")){
      setwd("./exprSet_ann")
    } else {
      dir.create("exprSet_ann")
      setwd("./exprSet_ann")    
    }
    if(dir.exists(gse)){
      setwd(gse)
    } else {
      dir.create(gse)
      setwd(gse)    
    }
    write_tsv(exprSet, paste(gse, "_matrix_expr.txt",sep=""))
    write_tsv(final_exprSet, paste(gse, "_matrix_expr_symbol.txt",sep=""))
    file.copy(paste(basedir,gse,fileToName,sep="/"), "./")
    write_tsv(pData(matrixData[[1]]), paste(gse, "pdata.txt",sep="_"))
    if(length(which(is.na(final_exprSet[-1])))>0){
      print(paste("There is NA in", gse, "using knn"))
      final_exprSet_noNA <- cbind(final_exprSet[1],as.data.frame(impute.knn(as.matrix(final_exprSet[-1]))$data))
      write_tsv(final_exprSet_noNA, paste(gse, "_matrix_expr_symbol_noNA.txt",sep = ""))
      # no dup
      final_exprSet_noNA_nodup <- final_exprSet_noNA %>% group_by_at(1) %>% summarise_all(funs(mean))
      write_tsv(final_exprSet_noNA_nodup, paste(gse, "_matrix_expr_symbol_noNAnoDup.txt",sep = ""))
    } else {
      final_exprSet_noNA_nodup <- final_exprSet %>% group_by_at(1) %>% summarise_all(funs(mean))
      write_tsv(final_exprSet_noNA_nodup, paste(gse, "_matrix_expr_symbol_noNAnoDup.txt",sep = ""))
    }
    dirNow <- getwd()
    setwd(basedir)
    write_lines(paste(dirNow,"/",gse, "_matrix_expr_symbol_noNAnoDup.txt",sep = ""), "../result.txt", append = TRUE)
  }
}