# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we add cell-type flags to denote where a given cell type may be located

library(Seurat)
library(tidyverse)

# Import path pointers

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples)

# First for each slide we will create metavariables that flag the location of a cell-type in a spot

c2l_assay <- "c2l_props"

# Eveything except cardiomyocytes and fibroblast must represent 10% of celltype score
ct_prop_param <- tibble(cts = c("cardiomyocyte", "fibroblast", "endothelial")) %>%
  mutate(prop_param = 0.125)

add_ct_flags <- function(slide, ct_prop_param, cell_props) {
  
  for(ct in ct_prop_param$cts) {
    
    ix <- grepl(ct, ct_prop_param$cts)
    
    prop_param <- ct_prop_param$prop_param[ix]
    
    slide[[paste0(ct,"_flag")]] = ifelse(cell_props[, ct] <= prop_param, 0, 1)
    
  }
  
  return(slide)
  
}

walk(visium_df$visium_file, function(visium_file) {
  
  print(visium_file)
  
  slide <- readRDS(visium_file)
  cell_props <- GetAssayData(slide, assay = c2l_assay) %>% t()
  
  slide_ct_prop_param <- ct_prop_param %>%
    dplyr::filter(cts %in% colnames(cell_props))
  
  slide <- add_ct_flags(slide = slide,
                        ct_prop_param = slide_ct_prop_param,
                        cell_props = cell_props)
  
  saveRDS(slide, file = visium_file)
  
})



