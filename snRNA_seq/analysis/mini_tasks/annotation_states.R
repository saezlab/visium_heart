# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we manually create annotation dictionaries
#' 
#' After integrating data, generating major annotations and
#' putative states, we will first create a unique dictionary of cell_type annotations
#' to create an object ready to be used for deconvolution via SPOTlight or
#' for the conversion mirror for cell2location
#' 

library(Seurat)
library(tidyverse)

deconv_col <- "deconv_state"

# Putative CM states
scell_data <- readRDS("./visium_results_manuscript/ct_data/cardio_snRNA.Rds")
vars_to_transfer <- "States"
scell_data_meta <- scell_data@meta.data %>%
  rownames_to_column("cell_id") %>% 
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id, split = "-"), ~.x[1])) %>%
  dplyr::select(orig.ident, cell_id, all_of(vars_to_transfer)) %>%
  dplyr::rename(deconv_col = vars_to_transfer)
  
scell_data_meta_cm <- scell_data_meta

# Putative fibro states
scell_data <- readRDS("./visium_results_manuscript/ct_data/fibro_coembd/fibroblasts_states.rds")
vars_to_transfer <- "States"
scell_data_meta <- scell_data@meta.data %>%
  rownames_to_column("cell_id") %>% 
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id, split = "-"), ~.x[1])) %>%
  dplyr::select(orig.ident, cell_id, all_of(vars_to_transfer)) %>%
  dplyr::rename(deconv_col = vars_to_transfer)

scell_data_meta_fb <- scell_data_meta

rm(scell_data)

# All_ans
all_annotations <- rbind(scell_data_meta_cm, 
                         scell_data_meta_fb)

write.table(all_annotations, 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = "\t",
            file = "./visium_results_manuscript/ct_data/state_annotations.txt")

























