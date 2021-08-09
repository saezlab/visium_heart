# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we read Hubners atlas
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(optparse)
library(scater)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "path where the folder with a collection of Seurat objects is"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds objects")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Evaluate parameters -----------------------------------------------
slide <- readRDS(path)
print("reading worked")

meta_hca <- slide@meta.data

meta_hca <- meta_hca %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(region == "LV",
                Used == "Yes",
                !cell_type %in% c("NotAssigned","doublets", "nan"))

# Now we need 3 dsets nuclei, cells, nuclei + cells, and both cell_type and cell_states
nuclei_meta_hca <-  meta_hca %>%
  dplyr::filter(source %in% c("CD45+", "Nuclei"))

cell_meta_hca <-  meta_hca %>%
  dplyr::filter(source %in% c("Cells"))

# Filtering annotated useful cells ---------------------------------
all_slide <- slide@assays$RNA@counts[ , meta_hca$cell_id]

nuclei_slide <- slide@assays$RNA@counts[ , nuclei_meta_hca$cell_id]

cell_slide <- slide@assays$RNA@counts[ , cell_meta_hca$cell_id]

# Function to generate pseudobulk counts --------------------------

get_pseudobulk_info <- function(seurat_obj, meta) {
  
  vars <- set_names(c("cell_type", "cell_states"))
  
  # Creates pseudobulk profile for each var
  bulk_p_data <- map(vars, function(x) { 
    sumCountsAcrossCells(x = seurat_obj,
                         ids = meta[, x])
  })
  
  bulk_p_data[["annotations"]] <- meta
  
  return(bulk_p_data)
}

# Final object ----------------------------------------------------

pseudobulk_info <- list("all" = get_pseudobulk_info(all_slide, meta_hca),
                        "nuclei" = get_pseudobulk_info(nuclei_slide, nuclei_meta_hca),
                        "cells" = get_pseudobulk_info(cell_slide, cell_meta_hca))

saveRDS(pseudobulk_info, out_path)