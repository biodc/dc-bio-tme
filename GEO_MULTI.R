target_geo <- FALSE
target_geo_file <- read.table("./get_target.txt", header = TRUE)
target_tcga <- FALSE
do_ciber <- TRUE
tcga_project <- "TCGA-STAD"
tcga_experimental_strategy <- "RNA-Seq"
tcga_data_category <- "Transcriptome Profiling"
tcga_data_type <- "Gene Expression Quantification"
tcga_workflow_type <- "HTSeq - FPKM"
if(!dir.exists("workspace")){
  dir.create("./workspace")
}
if(!file.exists("./result.txt")){
  file.remove("./result.txt")
}
if(target_geo) {
  setwd("./workspace")
  cores <- detectCores()/2
  c <- makeCluster(cores)
  parApply(c, target_geo_file,1, function(r){
    source("../Required_packages.R")
    gse <- r[[1]]
    t <- r[[2]]
    if(t=="A"){
      source("../affy_process.R")
      process_geo_affy(gse)
    } else if (t=="I") {
      source("../beadarray_process.R")
      process_geo_ill(gse)
    }
    return("Success!")
  })
  stopCluster(c)
  setwd("..")
}
if(target_tcga) {
  setwd("./workspace")
  source("../Required_packages.R")
  source("../get_tcga.R")
  query <- GDCquery(project = tcga_project, 
                    data.category = tcga_data_category, 
                    data.type = tcga_data_type, 
                    workflow.type = tcga_workflow_type, 
                    experimental.strategy = tcga_experimental_strategy)
  GDCdownload(query)
  tcgaData <- GDCprepare(query)
  process_tcga(tcgaData)
  setwd("..")
}
if(do_ciber){
  source("./Required_packages.R")
  fileList <- read_lines("./result.txt")
  if(!dir.exists("ciber")){
    dir.create("ciber")
  }
  setwd("ciber")
  sapply(fileList, function(r){
    source("../process_ciber.R")
    do_cibersort(r)
    return("Success!")
  })
  setwd("..")
}

