# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform differential expression analysis

library(optparse)
library(tidyverse)
library(Seurat)
source("./utils/dea.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--out_df"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds object with differential features"),
  make_option(c("--group_class"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "Identity to characterize"),
  make_option(c("--test_assays"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "assays to test separated by commas"),
  make_option(c("--lfc"), 
              action= "store", 
              default = "0.05", 
              type = 'character',
              help = "log fold change threshold"),
  make_option(c("--only_pos"), 
              action= "store", 
              default = "yes", 
              type = 'character',
              help = "log fold change threshold")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Read object to analyze and define assay -----------------------------------------------------
integrated_data <- readRDS(data_path)
Idents(integrated_data) <- group_class

# Define_assays -----------------------------------------------------
test_assays <- unlist(strsplit(test_assays, ","))

# Positive?

if(only_pos == "yes") { 
  only_pos = TRUE
} else {
  only_pos = FALSE
}

# Run differential test
dea_res <- find_allfeat(visium_slide = integrated_data, 
             assays_collection = test_assays,
             logfc.threshold = as.numeric(lfc),
             only.pos = only_pos)

saveRDS(dea_res, 
        file = out_df)
