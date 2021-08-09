# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we create a dictionary of gene_sets associated with states
#' The differeentially expressed genes are coming from the
#' funcomics pipeline

library(Seurat)
library(tidyverse)

dea_folder <- "./visium_results_manuscript/ct_data/dea/"

dea_files <- list.files("./visium_results_manuscript/ct_data/dea",
                        full.names = F)
cts <- map_chr(strsplit(dea_files, "_dea"), ~.x[1])

param_df <- tibble(dea_file = paste0(dea_folder, dea_files), 
                   cell_type = cts) %>%
  mutate(dea_obj = map(dea_file, readRDS)) %>%
  mutate(dea_obj = map(dea_obj, ~ .x[["RNA"]] %>%
                         group_by(cluster) %>%
                         dplyr::filter(p_val_adj < 0.00001,
                                       avg_logFC >= 0.5))) %>%
  unnest() %>%
  dplyr::select(cell_type, p_val_adj, 
                avg_logFC, gene, cluster) %>%
  dplyr::mutate(cell_state = paste0(cell_type, "_", cluster)) %>%
  arrange(cell_state, -avg_logFC, gene) %>%
  dplyr::select(cell_state, gene) %>%
  group_by(cell_state)
  
gsets <- param_df %>% 
  nest() %>% 
  mutate(data = map(data, ~.x[[1]])) %>% 
  deframe() %>%
  saveRDS(file = "./visium_results_manuscript/ct_data/dea/cellstate_gsets.rds")























