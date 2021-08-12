# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' For each slide, we add module scores of Niches to 
#' be modelled in MISTy
#' 

library(tidyverse)
library(Seurat)
source("./analysis/utils/funcomics.R")

visium_folder <- "./processed_visium/objects/"

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples)

niche_markers <- readRDS("./markers/niche_markers.rds")

# Add niche scores --------------------------------
add_niche_assay <- function(visium_file, sample) {
  
  # read visium and table
  print(sample)
  
  visium_slide <- readRDS(visium_file)
  
  visium_slide <- getTF_matrix_MS(visium_slide,
                                  MS_regulon = niche_markers,
                                  assay = "SCT",
                                  module_name = "niche_footprint")
  
  saveRDS(visium_slide, file = visium_file)
  
  DefaultAssay(visium_slide) = "niche_footprint"
  
  SpatialFeaturePlot(visium_slide, features = rownames(visium_slide))
}

pwalk(visium_df, add_niche_name)











