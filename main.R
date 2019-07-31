all_green <- TRUE

target_geo <- FALSE
target_tcga <- FALSE
merge_gene_true <- FALSE

batch_correction <- FALSE
do_ciber <- FALSE

purge_result <- FALSE
purge_workspace <- FALSE

do_mapCounter <- FALSE
do_merge_map_ciber <- FALSE

do_ciber_map_cluster <- TRUE

tcga_project <- "TCGA-STAD"
tcga_experimental_strategy <- "RNA-Seq"
tcga_data_category <- "Transcriptome Profiling"
tcga_data_type <- "Gene Expression Quantification"
tcga_workflow_type <- "HTSeq - FPKM"

if (all_green) {
  target_geo <- TRUE
  target_tcga <- TRUE
  merge_gene_true <- TRUE
  batch_correction <- TRUE
  do_ciber <- TRUE
  purge_result <- TRUE
  purge_workspace <- TRUE
  do_mapCounter <- TRUE
  do_merge_map_ciber <- TRUE
  do_ciber_map_cluster <- TRUE
}
ciber_matrix <- "./output/all_gene_matrix_after_batch_correction.txt"

ciber_result <- "./output/cibersort_all_result.txt"
mcp_result <- "./output/mcp_all_result.txt"

source("./func/init_workspace.R")
init_workspace()

if (target_geo) {
  source("./func/process_geo.R")
  process_geo()
}

if (target_tcga) {
  source("./func/process_tcga.R")
  process_tcga()
}

if (merge_gene_true) {
  source("./func/merge_gene_matrix.R")
  merged_gene_matrix <- mergeGeneMatrix("./result.txt")
}

if (batch_correction) {
  source("./func/batch_correction.R")
  mabc <- process_batch_correction()
}

if (do_ciber) {
  source("./func/process_ciber.R")
  if (exists("mabc")) {
    process_ciber(mabc)
  } else {
    process_ciber(ciber_matrix)
  }
}

if (do_mapCounter) {
  source("./func/process_mcp.R")
  if (exists("mabc")) {
    process_mcp(mabc)
  } else {
    process_mcp(ciber_matrix)
  }
}

if (do_ciber_map_cluster) {
  source("./func/process_merge_map_ciber.R")
  merge_map_ciber(mcp_result, ciber_result)
}



