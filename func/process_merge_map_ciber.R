merge_map_ciber <- function(mcpPath, ciberPath) {
  ciberData <- read_tsv(ciberPath)
  mcpData <- read_tsv(mcpPath)
  ciberData %<>% column_to_rownames("immCell") %>% t %>% as.data.frame %>% rownames_to_column("Sample") %>% filter(`P-value`<0.05) %>% column_to_rownames("Sample") %>% select_at(1:22) %>% rownames_to_column("Sample")
  mcpData %<>% column_to_rownames("Type") %>% t %>% as.data.frame %>% select_at(9) %>% rownames_to_column("Sample")
  final_matrix <- inner_join(ciberData, mcpData, by="Sample") %>% column_to_rownames("Sample") %>% t
  save(final_matrix, file="./output/map_ciber_final_matrix.rda")
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
  setwd("..")
}