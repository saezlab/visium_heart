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
              help = "path pointing to the folder where your sce object will be created")
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
sce_file <- paste0(sce_file_root, file_name)

cell_data <- readRDS(seurat_file)
cell_data <- as.SingleCellExperiment(cell_data)

saveHDF5SummarizedExperiment(cell_data, dir = sce_file, prefix = "", replace = FALSE,
                             chunkdim = NULL, level = NULL, as.sparse = NA,
                             verbose = NA)



