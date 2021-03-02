# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we transform data from scell objects to anndatas

library(Seurat)
library(zellkonverter)
library(tidyverse)
library(optparse)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--scell_data"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell data with states in a variable"),
  make_option(c("--out_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the anndata object?")
  
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

scell_obj <- readRDS(scell_data)
# Replace all dots with _
colnames(scell_obj@meta.data) <- gsub("[.]", "_", colnames(scell_obj@meta.data))
scell_obj <- as.SingleCellExperiment(scell_obj, assay = "RNA")
writeH5AD(scell_obj, file = out_file)










