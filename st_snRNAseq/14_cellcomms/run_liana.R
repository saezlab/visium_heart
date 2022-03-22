# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we select from a single integrated object
#' cell-types of interest. Usually a single cell is expected
#' 
#' Once a selection is done, a whole integration process
#' and optimization of clustering is done to
#' find cell states.
#' 

library(optparse)
library(tidyverse)
library(Seurat)
library(liana)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--out_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the liana results"),
  make_option(c("--cell_class"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "cell classes to consider (separated by commas)")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Cell types to include ----------------------------------------------------------------------------
cell_type <- unlist(strsplit(cell_class, ","))

# Read object to subset -----------------------------------------------------
integrated_data <- readRDS(data_path)

DefaultAssay(integrated_data) <- "RNA"

true_ix <- integrated_data@meta.data[, "cell_type"] %in% cell_type

# subset based on filtering and quickly get characteristic profile
integrated_data <- integrated_data[ , true_ix]

# Set identities
Idents(integrated_data) <- "annotation"

# Run liana
liana_test <- liana_wrap(integrated_data, assay = "RNA")

# Save results
saveRDS(liana_test, file = out_file)

# Save object
saveRDS(integrated_data, file = gsub(pattern = "[.]rds",
                                     replacement = "_obj.rds", 
                                     x = out_file))











