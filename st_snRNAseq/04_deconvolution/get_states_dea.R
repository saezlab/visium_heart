# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we create a dictionary of gene_sets associated with states
#' The differeentially expressed genes are coming from the
#' funcomics pipeline

library(tidyverse)

# These are the markers of the cell-states
state_mrkrs <- tibble(marker_file = list.files("./cell_states", full.names = T)) %>%
  dplyr::mutate(cell_type = gsub("./cell_states/", "", marker_file)) %>%
  dplyr::mutate(marker_file = paste0(marker_file,"/annotation.rds")) %>%
  dplyr::mutate(markers = map(marker_file, readRDS)) %>%
  dplyr::select(cell_type, markers) %>%
  unnest() %>%
  dplyr::select(cell_type, p_val_adj, cluster, gene, avg_log2FC) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(-p_val_adj) %>%
  mutate(source = paste0(cell_type, "_", cluster),
         mor = sign(avg_log2FC)) %>%
  dplyr::rename("target" = gene,
                "likelihood" = avg_log2FC) %>%
  dplyr::select(-c("cell_type", "cluster"))

# I will ignore all state mrkrs that have less than 50 markers
state_mrkrs <- state_mrkrs %>%
  group_by(source) %>%
  nest() %>%
  dplyr::mutate(n_mrkrs = map(data, nrow)) %>%
  #dplyr::filter(n_mrkrs >= 50) %>%
  dplyr::select(-n_mrkrs) %>%
  dplyr::filter(! grepl("Adipo", source),
                ! grepl("Mast", source)) %>%
  unnest() %>%
  ungroup()

saveRDS(state_mrkrs, file = "./results/ct_data/state_genesets.rds")

# Make a list of markers for the module score estimation

state_mrkrs_list <- state_mrkrs %>%
  dplyr::select(source, target) %>%
  dplyr::group_by(source) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

saveRDS(state_mrkrs_list, file = "./results/ct_data/state_genesets_list.rds")

