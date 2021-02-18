# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we manually annotate clusters to be used in subsequent
#' analyses
#' 
#' Warning: Hard coding is expected here, run plot_knownmarkers.R after this
#' for validation
#' 

library(optparse)
library(tidyverse)
library(Seurat)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--dictionary_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "tsv file with annotations to transfer (only one extra column)"),
  make_option(c("--object_id"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "left id for left join"),
  make_option(c("--dictionary_id"),
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "right id for left join"),
  make_option(c("--new_variable"),
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "new variable in meta.data from dictionary_id"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds object")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Read object to annotate and get meta-data
scell_obj <- readRDS(data_path)
meta_data <- scell_obj@meta.data
# I assume it will always be a numeric cluster label
meta_data[, object_id] <- as.numeric(as.character(meta_data[,object_id] ))

# Read object to annotate and get meta-data
annotation_data <- read.table(file = dictionary_path,
                              header = T,
                              sep = "\t",
                              stringsAsFactors = F)

# Left-join
merge_conditional <- set_names(dictionary_id,object_id)

meta_data <- left_join(meta_data,
                       annotation_data,
                       by = merge_conditional)

# Add info to object
scell_obj <- AddMetaData(scell_obj, 
                         metadata = meta_data[, new_variable],
                         col.name = new_variable)

# Save
saveRDS(scell_obj, file = out_path)



































































































