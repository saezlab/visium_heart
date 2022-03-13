# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we select from a single integrated object
#' cell-types of interest. Usually a single cell is expected
#' 
#' Once a selection is done, liana is run using the annotation column
#' we can overwrite annotations to keep regular cell-types if needed
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
  make_option(c("--cell_type_list"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "cell types to consider (separated by commas)"),
  make_option(c("--cell_state_list"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "cell states to consider (separated by commas)")
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
ct_list <- unlist(strsplit(cell_type_list, ","))
cs_list <- unlist(strsplit(cell_state_list, ","))

cell_type <- c(ct_list, cs_list)

# Read object to subset -----------------------------------------------------
integrated_data <- readRDS(data_path)

DefaultAssay(integrated_data) <- "RNA"

true_ix <- integrated_data@meta.data[, "cell_type"] %in% cell_type

# subset based on filtering and quickly get characteristic profile
integrated_data <- integrated_data[ , true_ix]

# make the trick so that you can analyze cell-types and cell-states
copy_meta <- integrated_data@meta.data %>%
  as.data.frame() %>%
  mutate(annotation = as.character(annotation)) %>%
  mutate(annotation = ifelse(cell_type %in% ct_list, cell_type, annotation))

integrated_data$annotation <- as.factor(as.character(copy_meta$annotation))

# Set identities
Idents(integrated_data) <- "annotation"

# Run liana
liana_test <- liana_wrap(integrated_data, assay = "RNA")

# Save results
saveRDS(liana_test, file = out_file)











