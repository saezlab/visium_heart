# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Convert from Seurat to anndata
#' 
#' 

library(tidyverse)
library(Seurat)
library(loomR)
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

DefaultAssay(scell_obj) <- "RNA"

scell_obj_loom <- as.loom(scell_obj, 
                          filename = out_file, 
                          verbose = FALSE)


