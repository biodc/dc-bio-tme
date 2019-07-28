merge_map_ciber <- function(mcpPath, ciberPath) {
  ciberData <- read_tsv(ciberPath)
  mcpData <- read_tsv(mcpPath)
  trim_ciberData <- (column_to_rownames(ciberData, "immCell"))
  trim_ciberData <- trim_ciberData[!rownames(trim_ciberData) %in% c("P-value", "Correlation", "RMSE"),]
  t_trim_ciberData <- rownames_to_column(as.data.frame(t(trim_ciberData)), "Sample")
  trim_mcpData <- (column_to_rownames(mcpData, "Type"))
  trim_mcpData <- trim_mcpData[rownames(trim_mcpData)=="Fibroblasts",]
  t_trim_mcpData <- rownames_to_column(as.data.frame(t(trim_mcpData)), "Sample")
  final_matrix <- t(column_to_rownames(full_join(t_trim_ciberData, t_trim_mcpData, by="Sample"), "Sample"))
  save(final_matrix, file="map_ciber_final_matrix.rda")
  final_matrix_to_save <- rownames_to_column(as.data.frame(final_matrix), "Type")
  write_tsv(final_matrix_to_save, "./output/map_ciber_final_matrix.txt")
  final_matrix_norm <- preprocessCore::normalize.quantiles(final_matrix)
  rownames(final_matrix_norm) <- rownames(final_matrix)
  colnames(final_matrix_norm) <- colnames(final_matrix)
  setwd("./output/")
  ConsensusClusterPlus::ConsensusClusterPlus(d = final_matrix_norm, 
                                             maxK = 5, 
                                             reps = 1000, 
                                             pItem = 0.8,
                                             pFeature = 1,
                                             title = "ClusterResult",
                                             innerLinkage = "ward.D2",
                                             finalLinkage = "ward.D2",
                                             clusterAlg = "hc",
                                             distance = "euclidean",
                                             writeTable = TRUE,
                                             plot="png",
                                             )
}