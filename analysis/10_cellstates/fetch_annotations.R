# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Get joint annotation of cell-states

library(Seurat)
library(tidyverse)

folder = "./atac_rna_states"
all_objs <- list.files(folder, recursive = T,full.names = T)

all_annotations <- map(all_objs, function(fname) {
  
  cell_obj <- readRDS(fname)
  
  if("tech" %in% colnames(cell_obj@meta.data)) {
    
    cell_obj@meta.data %>% 
      as.data.frame() %>%
      rownames_to_column("raw_id") %>%
      dplyr::filter(tech == "RNA") %>%
      dplyr::select(raw_id, orig.ident,cell_type, annotation)
    
  } else {
    
    cell_obj@meta.data %>% 
      as.data.frame() %>%
      rownames_to_column("raw_id") %>%
      dplyr::select(raw_id, orig.ident,cell_type, annotation)
    
  }
  
}) 

saveRDS(all_annotations, "./processed_snrnaseq/cell_states/cellstate_annotation_list.rds")
    




