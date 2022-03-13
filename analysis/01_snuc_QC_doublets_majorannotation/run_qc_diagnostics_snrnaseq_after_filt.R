# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we recover the QC stats of all slides
#' after individual processing

library(tidyverse)
library(Seurat)

folder = T
path = "./processed_snrnaseq/objects/"
sample_names <- list.files(path)
sample_names <- gsub("[.]rds", "", sample_names)
slide_files <- set_names(paste0(path, sample_names, ".rds"), sample_names)

qc_afetr_filt <- map(slide_files, function(rds_file) {
  print(rds_file)
  visium_slide_meta <- readRDS(rds_file)@meta.data %>%
    group_by(orig.ident) %>%
    summarise(n_cell_after_filt = length(orig.ident),
              median_counts = median(nCount_RNA),
              median_ngenes = median(nFeature_RNA))
  
  return(visium_slide_meta)
  
}) %>%
  enframe("sample_id") %>%
  unnest() %>%
  dplyr::select(-orig.ident) %>%
  write.table(row.names = F, col.names = T, quote = F, sep = ",",
              file = "./processed_snrnaseq/initial_qc/all_qcs_after_processing.csv")
