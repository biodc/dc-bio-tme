init_workspace <- function() {
  source("./func/Required_packages.R")
  if(!dir.exists("workspace")){
    dir.create("./workspace")
  }
  if(!dir.exists("output")){
    dir.create("./output")
  }
  if(purge_result && file.exists("./result.txt")){
    file.remove("./result.txt")
  }
  if(purge_workspace && dir.exists("./workspace")) {
    unlink("workspace", recursive = TRUE)
    dir.create("workspace")
  } 
}