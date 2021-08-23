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
  make_option(c("--file_name"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "name of the seurat object without extension (we assume is .rds)"),
  make_option(c("--seurat_file_root"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "path pointing to the folder where your seurat object is located"),
  make_option(c("--sce_file_root"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "path pointing to the folder where your sce object will be created"),
  make_option(c("--subsample"), 
              action = "store_true", 
              default = FALSE, 
              type = 'logical',
              help = "subset the data file?"),
  make_option(c("--subsample_group"), 
              action ="store", 
              default = "orig.ident", 
              type = 'character',
              help = "variable used to group objects"),
  make_option(c("--subsample_prop"), 
              action ="store", 
              default = 0.5, 
              type = 'double',
              help = "proportion per grouping variable to keep")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

seurat_file <- paste0(seurat_file_root, file_name, ".rds")
sce_file <- paste0(sce_file_root, file_name, ".rds")

set.seed(241099)

if(subsample) {
  
  cell_data <- readRDS(seurat_file)
  #First we need to count how many samples per group we are gonna pick
  
} else {
  
  cell_data <- readRDS(seurat_file)
  
  cell_data <- DietSeurat(
    cell_data,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = "RNA",
    dimreducs = "umap"
  )
  
  cell_data <- cell_data[1:100, 1:100]
  
  cell_data <- as.SingleCellExperiment(cell_data)
  
}

saveRDS(object = cell_data,
        compress = F,
        file = sce_file)

#saveHDF5SummarizedExperiment(cell_data, dir = sce_file, prefix = "", replace = FALSE,
#                             chunkdim = NULL, level = NULL, as.sparse = NA,
#                             verbose = NA)



