# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we generate the state and cell-type matrices to characterize niches

library(tidyverse)
library(Seurat)
source("./analysis/utils/misty_pipeline.R")

# First cell-type proportions

# Main ------------------------------------------------------------------------
# Getting sample annotations --------------------------------------------------
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# First cell-types

param_df <- tibble(slide = slide_ids, 
                   slide_file = paste0(slide_files_folder, slide_files)) %>%
  dplyr::mutate(state_assay = map2(slide, slide_file, function(s, sf) {
    
    print(s)
    
    vs <- readRDS(sf)
    
    s_ass <- GetAssayData(vs, assay = "c2l_props")
    
    colnames(s_ass) <- paste0(s, "..", colnames(s_ass))
    
    s_ass
    
  }))

# Put all compositions in a single object

cellprops_mat <- reduce(param_df$state_assay, cbind)

cellprops_mat <- t(cellprops_mat)

cellprops_info <- cellprops_mat %>%
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  pivot_longer(-spot_id)

saveRDS(cellprops_info, file = "./results/ind_mats/cell_type_props.rds")

# Second cell-states

param_df <- tibble(slide = slide_ids, 
                   slide_file = paste0(slide_files_folder, slide_files)) %>%
  dplyr::mutate(state_assay = map2(slide, slide_file, function(s, sf) {
    
    print(s)
    
    vs <- readRDS(sf) %>%
      positive_states(., assay = "cell_states") %>%
      filter_states(slide = .,
                    by_prop = F,
                    prop_thrsh = 0.1)
    
    s_ass <- GetAssayData(vs, assay = "cell_states_pos")
    
    colnames(s_ass) <- paste0(s, "..", colnames(s_ass))
    
    s_ass
    
  }))

# Put all states in a single object

state_score_mat <- reduce(param_df$state_assay, cbind)

state_score_mat <- t(state_score_mat)

rm(param_df)

# Per state
state_score_mat <- state_score_mat %>%
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  pivot_longer(-spot_id) %>%
  dplyr::filter(!grepl("Adipo", name)) %>%
  dplyr::filter(!grepl("Mast", name))

saveRDS(state_score_mat, file = "./results/ind_mats/cell_state_scorespos.rds")



