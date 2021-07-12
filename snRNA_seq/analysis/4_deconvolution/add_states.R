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
gsets <- readRDS(file = "./visium_results_manuscript/ct_data/dea/cellstate_gsets.rds")

walk(visium_df$visium_file, function(visium_file) { 
  print(visium_file)
  
  visium_slide <- readRDS(visium_file) %>%
    getTF_matrix_MS(visium_slide = .,
                    MS_regulon = gsets,
                    assay = "SCT",
                    module_name = "cell_states")
  
  saveRDS(visium_slide, file = visium_file)
})

