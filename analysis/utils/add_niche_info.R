# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Integrate niche annotations to individual slides for simplification
#' 
#' visium_folder
#' |
#' -------samplename.rds
#' 
#' 
#' cpseudobulk from integrated slides
#' 

library(optparse)
library(tidyverse)
library(Seurat)
library(cluster)
library(cowplot)
source("./analysis/utils/funcomics.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--visium_folder"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "folder with visium slides"),
  make_option(c("--pseudobulk_file"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "name of the sample if single path provided")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Get integrated annotations
integrated_meta <- readRDS(pseudobulk_file)[[1]][["annotations"]] %>%
  rownames_to_column("spot_id") %>% 
  dplyr::mutate(spot_id = map_chr(strsplit(spot_id, "-"), ~ .x[[1]])) %>%
  dplyr::select(spot_id, orig.ident, opt_clust_integrated) %>%
  dplyr::mutate_at(c("spot_id", "orig.ident", "opt_clust_integrated"),
                   as.character()) %>%
  dplyr::mutate(opt_clust_integrated = paste0("niche_", opt_clust_integrated))

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  dplyr::filter(sample %in% pull(integrated_meta, 
                                 orig.ident) %>% unique())

add_niche_name <- function(visium_file, sample) {
  
  # read visium and table
  print(sample)
  
  visium_slide <- readRDS(visium_file)
  
  sample_meta <- integrated_meta %>%
    dplyr::filter(orig.ident == sample)
  
  visium_meta <- visium_slide@meta.data %>%
    rownames_to_column("spot_id") %>% 
    dplyr::mutate(spot_id = map_chr(strsplit(spot_id, "-"), ~ .x[[1]])) %>%
    mutate_at(c("spot_id"), as.character()) %>%
    left_join(sample_meta, by = "spot_id")
  
  visium_slide$opt_clust_integrated <- visium_meta$opt_clust_integrated
    
  saveRDS(visium_slide, file = visium_file)
}

pmap(visium_df, add_niche_name)
