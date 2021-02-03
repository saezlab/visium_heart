# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create pseudobulk profiles from a directory with
#' Seurat data sets.
#' The input of this program requires to define for each data set a column
#' with identity

library(tidyverse)
library(Seurat)
library(optparse)
library(scater)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "path where the folder with a collection of Seurat objects is"),
  make_option(c("--vars"), 
              action = "store", 
              default = NULL, 
              type = 'character',
              help = "what variables to use as groups"),
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
slide_files <- list.files(path)
# Assume rds files only
slide_files <- slide_files[grepl("[.]rds", slide_files)]
slide_files <- paste0(path, slide_files)

# Identify levels of pseudobulking

vars <- set_names(unlist(strsplit(vars, ",")))

# Read data and create pseudobulk profiles at two levels ----------------------------------

get_sample_pseudo <- function(slide_file, vars) {
  
  slide <- readRDS(slide_file)
  
  # This is only to deal with the lack of annotations in the meta_data
  slide$cell_type <- gsub(" ", "_", (Idents(slide)))
  slide$major_cell_type <- gsub(" ", "_", slide$top_annotation) 
  
  # Creates pseudobulk profile for each var
  bulk_p_data <- map(vars, function(x) { 
    sumCountsAcrossCells(x = as.matrix(slide@assays$RNA@counts),
                         ids = slide@meta.data[, x])
  })
  
  bulk_p_data[["annotations"]] <- unique(slide@meta.data[,vars])
  
  return(bulk_p_data)
}

pseudobulk_profiles <- map(set_names(slide_files), 
                           get_sample_pseudo, 
                           vars = vars)

saveRDS(pseudobulk_profiles, 
        file = out_path)











































































