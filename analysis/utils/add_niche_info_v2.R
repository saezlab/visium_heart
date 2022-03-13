# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Integrate niche annotations to individual slides for simplification
#' 
#' Coming from the cell-type based niche definition
#' 
#' visium_folder
#' |
#' -------samplename.rds
#' 
#' 
#' annotation table from find_niches_ct
#' 

library(Seurat)
library(tidyverse)

spot_annotation <- readRDS("./results/niche_mapping/ct_niches/niche_annotation_ct.rds") %>%
  dplyr::mutate(spot_id = map_chr(strsplit(row_id, "[..]"), ~.x[[3]]),
                orig.ident = map_chr(strsplit(row_id, "[..]"), ~.x[[1]])) %>%
  dplyr::select_at(c("orig.ident", "spot_id", "ct_niche")) %>%
  dplyr::rename("opt_clust_integrated_ct" = ct_niche)

# Get visium slides --------------------------------
visium_folder <- "./processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  dplyr::filter(sample %in% pull(spot_annotation, 
                                 orig.ident) %>% unique())

# Add the niche information -------------------------

add_niche_name <- function(visium_file, sample) {
  
  # read visium and table
  print(sample)
  
  visium_slide <- readRDS(visium_file)
  
  sample_meta <- spot_annotation %>%
    dplyr::filter(orig.ident == sample)
  
  if(ncol(visium_slide) == nrow(sample_meta)) {
  
  visium_meta <- visium_slide@meta.data %>%
    rownames_to_column("spot_id") %>% 
    dplyr::select(spot_id, orig.ident) %>%
    mutate_at(c("spot_id"), as.character()) %>%
    left_join(sample_meta, by = "spot_id")
  
  # Here we overwrite
  visium_slide$opt_clust_integrated <- visium_meta$opt_clust_integrated_ct
  
  saveRDS(visium_slide, file = visium_file)
  
  } else{
    
    print("something went wrong with the mapping")
    
  }
}

pwalk(visium_df, add_niche_name)
