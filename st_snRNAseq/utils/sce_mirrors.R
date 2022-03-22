# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform the transformation needed for the 
#' single cell objects to generate the shiny apps
#' powered by iSEE

library(optparse)
library(Seurat)
library(HDF5Array)
library(scater)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--seurat_file"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell data with states in a variable"),
  make_option(c("--sce_folder"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the anndata object?"),
  make_option(c("--assay"), 
              action= "store", 
              default = "RNA", 
              type = 'character',
              help = "where to save the anndata object?"),
  make_option(c("--reduction"), 
              action= "store", 
              default = "umap", 
              type = 'character',
              help = "where to save the anndata object?")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

cell_data <- readRDS(seurat_file)
  
cell_data <- DietSeurat(
    cell_data,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = reduction
  )
  
cell_data <- as.SingleCellExperiment(cell_data)

saveHDF5SummarizedExperiment(cell_data, dir = sce_folder, 
                             prefix = "", replace = FALSE,
                             chunkdim = NULL, level = NULL, as.sparse = NA,
                             verbose = NA)
