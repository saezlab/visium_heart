# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Integrate cell2location outputs into single slide results
#' The expected format of the input folders looks like this:
#' 
#' visium_folder
#' |
#' -------samplename.rds
#' 
#' 
#' c2l_folder
#' |
#' -------samplename_W_cell_density_q05.csv
#' -------samplename_W_cell_density.csv
#' 

library(optparse)
library(tidyverse)
library(Seurat)
library(cluster)
library(cowplot)
source("./utils/funcomics.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--visium_folder"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "folder with visium slides"),
  make_option(c("--c2l_folder"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "name of the sample if single path provided"),
  make_option(c("--assay_name"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "name of the sample if single path provided"),
  make_option(c("--q05"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "use q05 matrix"),
  make_option(c("--make_proportions"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "use q05 matrix")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples)

# Get cell2location files --------------------------------
c2l_files <- list.files(c2l_folder, full.names = F)
if(q05) {
  c2l_files <- c2l_files[grepl("q05", c2l_files)]
} else {
  c2l_files <- c2l_files[!grepl("q05", c2l_files)]
}
c2l_samples <- map_chr(strsplit(c2l_files,"_"), 
                            ~ .x[1])

c2l_df <- tibble(c2l_file = paste0(c2l_folder, 
                                   c2l_files),
                    sample = c2l_samples) 

# main
param_df <- left_join(visium_df, c2l_df) %>%
  dplyr::select(sample, visium_file, c2l_file)


add_c2l_assay <- function(sample, visium_file, c2l_file) {
  
  # read visium and table
  print(sample)
  
  visium_slide <- readRDS(visium_file)
  mat <- read.csv(c2l_file,row.names = 1)
  rownames(mat) <- map_chr(strsplit(rownames(mat),"_"), ~ .x[2])
  
  if(q05) {
    colnames(mat) <- gsub("q05_spot_factors", "", colnames(mat))
  } else {
    colnames(mat) <- gsub("mean_spot_factors", "", colnames(mat))
  }
  
  # Ensure same spot order
  rownames(mat) <- map_chr(strsplit(rownames(mat),"-"), ~ .x[1])
  meta_visium <- visium_slide@meta.data %>%
    rownames_to_column("spot_id") %>%
    mutate(raw_spot_id = map_chr(strsplit(spot_id,"-"), ~ .x[1]))
  
  mat <- mat[meta_visium$raw_spot_id, ]
  rownames(mat) <- meta_visium$spot_id
  
  # adding assay
  visium_slide[[assay_name]] = CreateAssayObject(data = t(mat))
  
  # adding alternative assay
  if(make_proportions) {
    prop_mat <- base::apply(mat, 1, function(x) {
      
      x/sum(x)
      
    })
    
    visium_slide[[paste0(assay_name, "_props")]] = CreateAssayObject(data = prop_mat)
    
  }
  
  saveRDS(visium_slide, file = visium_file)
}

pmap(param_df, add_c2l_assay)
