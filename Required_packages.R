if (!requireNamespace("BiocManager", quietly = TRUE)){
  print("Bioconductor not installed, installing...")
  install.packages("BiocManager")
}
library(parallel)
if(!require(GEOquery)){
  print("GEOquery not installed, installing...")
  BiocManager::install("GEOquery")
  require(GEOquery)
}
if(!require(dplyr)){
  print("dplyr not installed, installing...")
  BiocManager::install("dplyr")
  require(dplyr)
}
if(!require(affy)) {
  print("affy not installed, installing...")
  BiocManager::install("affy")
  require(affy)
}
if(!require(limma)) {
  print("limma not installed, installing")
  BiocManager::install("limma")
  require(limma)
}
if(!require(annotate)){
  print("annotate not installed, installing")
  BiocManager::install("annotate")
  require(annotate)
}
if(!require(illuminaHumanv3.db)){
  print("illuminaHumanv3.db not installed, installing")
  BiocManager::install("illuminaHumanv3.db")
  require(illuminaHumanv3.db)
}
if(!require(org.Hs.eg.db)){
  print("org.Hs.eg.db not installed, installing")
  BiocManager::install("org.Hs.eg.db")
  require(org.Hs.eg.db)
}
if(!require(readr)){
  print("readr not installed, installing")
  install.packages("readr")
  require(readr)
}
if(!require(dplyr)){
  print("dplyr not installed, installing")
  install.packages("dplyr")
  require(dplyr)
}
if(!require(magrittr)){
  print("magrittr not installed, installing")
  install.packages(magrittr)
  require(magrittr)
}
if(!require(stringr)){
  print("stringr not installed, installing")
  install.packages("stringr")
  require(stringr)
}
if(!require(lumiHumanIDMapping)){
  print("lumiHumanIDMapping not installed, installing")
  BiocManager::install("lumiHumanIDMapping")
  require(lumiHumanIDMapping)
}
if(!require(impute)) {
  print("Impute not installed, installing...")
  BiocManager::install("impute")
  require("impute")
}