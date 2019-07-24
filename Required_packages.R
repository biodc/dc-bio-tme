if (!requireNamespace("BiocManager", quietly = TRUE)){
  print("Bioconductor not installed, installing...")
  install.packages("BiocManager")
}
require(parallel)
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
  BiocManager::install("readr")
  require(readr)
}
if(!require(magrittr)){
  print("magrittr not installed, installing")
  BiocManager::install(magrittr)
  require(magrittr)
}
if(!require(stringr)){
  print("stringr not installed, installing")
  BiocManager::install("stringr")
  require(stringr)
}
if(!require(impute)) {
  print("Impute not installed, installing...")
  BiocManager::install("impute")
  require("impute")
}
if(!require(e1071)){
  print("e1071 not installed, installing")
  install.packages("e1071")
  require(e1071)
}
if(!require(SummarizedExperiment)) {
  print("SummarizedExperiment not installed, installing...")
  BiocManager::install("SummarizedExperiment")
  require("SummarizedExperiment")
}
if(!require(preprocessCore)){
  print("preprocessCore not installed, installing")
  install.packages("preprocessCore")
  require(preprocessCore)
}
if(!require(TCGAbiolinks)) {
  print("TCGAbiolinks not installed, installing...")
  BiocManager::install("TCGAbiolinks")
  require("TCGAbiolinks")
}
