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

ct_path <- "./visium_results_manuscript/ct_data/"
scell_objs <- list.files("./visium_results_manuscript/ct_data/")
scell_objs <- scell_objs[grepl(".rds", scell_objs)]

input_tibble <- tibble(file = paste0(ct_path, scell_objs),
                       cell = gsub("_states[.]rds", "", scell_objs))

# Generate meta information tibble ------------------------------------
get_meta_data <- function(file, cell, vars_to_transfer = "opt_state") {
  scell_data <- readRDS(file)
  scell_data_meta <- scell_data@meta.data %>%
    rownames_to_column("cell_id") %>% 
    dplyr::mutate(cell_id = map_chr(strsplit(cell_id, split = "-"), ~.x[1])) %>%
    dplyr::select(orig.ident, cell_id, all_of(vars_to_transfer)) %>%
    dplyr::rename(deconv_col = vars_to_transfer) %>%
    dplyr::mutate(deconv_col = paste0(cell,"_",deconv_col))
}

input_tibble <- input_tibble %>%
  dplyr::mutate(meta_info = map2(file, cell, get_meta_data))

all_annotations <- reduce(input_tibble$meta_info,
                          bind_rows)


# In case VSMCs are not obtained ------------------------------------
# WARNING MANUAL annotation always is done
all_annotations <- all_annotations %>%
  dplyr::mutate(deconv_col = ifelse(deconv_col == "pericytes_1", 
                                    "VSMCs", 
                                    deconv_col))

# Write final table ------------------------------------
write.table(all_annotations, 
            col.names = T, 
            row.names = F,
            quote = F,
            sep = "\t",
            file = "./visium_results_manuscript/ct_data/state_annotations.txt")

# Putative CM states
scell_data <- readRDS("./visium_results_manuscript/ct_data/cardiomyocytes_states.rds")
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



























