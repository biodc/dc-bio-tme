target_geo <- read.table("./get_target.txt", header = TRUE)

cores <- detectCores()/2
c <- makeCluster(cores)
parApply(c, target_geo,1, function(r){
  source("./Required_packages.R")
  gse <- r[[1]]
  t <- r[[2]]
  if(t=="A"){
    source("./affy_process.R")
    process_geo_affy(gse)
  } else if (t=="I") {
    source("./beadarray_process.R")
    process_geo_ill(gse)
  }
  return("Success!")
})
stopCluster(c)