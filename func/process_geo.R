process_geo <- function(){
  target_geo_file <- read.table("./get_target.txt", header = TRUE)
  setwd("./workspace")
  cores <- detectCores()/2
  c <- makeCluster(cores)
  parApply(c, target_geo_file,1, function(r){
    source("../func/Required_packages.R")
    gse <- r[[1]]
    t <- r[[2]]
    if(t=="A"){
      source("../func/affy_process.R")
      process_geo_affy(gse)
    } else if (t=="I") {
      source("../func/beadarray_process.R")
      process_geo_ill(gse)
    }
    return("Success!")
  })
  stopCluster(c)
  setwd("..")
}