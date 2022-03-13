# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Add module score of gene sets of cellular states from
#' get_states_dea.R

library(Seurat)
library(tidyverse)
source("./analysis/utils/funcomics.R")

# Main

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()

# Add module score ---------------------------------------------
gsets <- readRDS(file = "./results/ct_data/state_genesets.rds")
gsets_list <- readRDS(file = "./results/ct_data/state_genesets_list.rds")

walk(visium_df$visium_file, function(visium_file) { 
  print(visium_file)
  
  visium_slide <- readRDS(visium_file) %>%
    get_wmean_score(visium_slide = .,
                    network = gsets,
                    assay = "SCT",
                    module_name = "cell_states")
  
  #visium_slide <- getTF_matrix_MS(visium_slide,
  #                                MS_regulon = gsets_list,
  #                                assay = "SCT",
  #                                module_name = "cell_states_ms")
  
  saveRDS(visium_slide, file = visium_file)
})

# Add GRN data

# cm_grn <- read_csv("./results/ct_regnets/CM/gene_cluster.csv") %>%
#   dplyr::select(gene, cluster) %>%
#   dplyr::rename("source" = cluster,
#                 "target" = gene) %>%
#   dplyr::mutate(likelihood = 1, 
#                 mor = 1,
#                 source = paste0("CM_", source))
# 
# walk(visium_df$visium_file, function(visium_file) { 
#   print(visium_file)
#   
#   visium_slide <- readRDS(visium_file) %>%
#     get_wmean_score(visium_slide = .,
#                     network = cm_grn,
#                     assay = "SCT",
#                     module_name = "cell_states_rn")
# 
#   saveRDS(visium_slide, file = visium_file)
# })
