# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we filter the differential expression analysis for annotation

library(tidyverse)
library(Seurat)

degs <- readRDS("./visium_results_manuscript/integration/integrated_data_degs.rds")

degs <- degs$RNA %>%
  dplyr::mutate(specificity = pct.1 - pct.2) %>%
  arrange(cluster, p_val_adj, -specificity) %>%
  dplyr::filter(avg_logFC > 0.5) %>%
  arrange(cluster, -specificity) %>%
  group_by(cluster) %>%
  dplyr::slice(1:20) %>%
  write.table(file = "./visium_results_manuscript/integration/integrated_data_mrkr_genes.txt",
              sep = "\t",
              row.names = F, 
              col.names = T,
              quote = F)

degs <- degs$RNA %>%
  dplyr::mutate(specificity = pct.1 - pct.2) %>%
  arrange(cluster, p_val_adj, -specificity) %>%
  dplyr::filter(avg_logFC > 0.5) %>%
  arrange(cluster, -specificity) %>%
  group_by(cluster) %>%
  dplyr::slice(1:20) %>%
  dplyr::filter(cluster == 7) %>%
  write.table(file = "./visium_results_manuscript/integration/integrated_data_mrkr_genes_clust7.txt",
              sep = "\t",
              row.names = F, 
              col.names = T,
              quote = F)









