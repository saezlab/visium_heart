# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we create a dictionary of gene_sets associated with states
#' The differeentially expressed genes are coming from the
#' funcomics pipeline

library(Seurat)
library(tidyverse)

dea_folder <- "./results/ct_data/"
ct_folders <- list.dirs(dea_folder,recursive = F,full.names = F)

# gene_list_suffix <- "state_genelist.rds"
# 
# state_df <- tibble(ct = ct_folders, gene_set_file = paste0(dea_folder, ct_folders, "/", gene_list_suffix)) %>%
#   mutate(gene_sets = map2(ct, gene_set_file, function(x,y) {
#     
#     set_list <- readRDS(y)
#     names(set_list) <- paste0(x, "_", names(set_list))
#     return(set_list)
#     
#   }))
#   
# saveRDS(state_df$gene_sets %>% unlist(recursive = F), file = "./results/ct_data/state_genesets.rds")
# 

param_df <- tibble(ct = ct_folders, 
       ct_folder = paste0(dea_folder, ct_folders, "/"))

gene_list_suffix <- "_state_mrks.txt"

get_de_genes <- function(ct_folder) {
  folder_files <- list.files(ct_folder)
  de_file <- folder_files[grepl(gene_list_suffix, folder_files)]
  de_table <- read_table2(paste0(ct_folder, de_file)) %>%
    dplyr::mutate(mt_genes = grepl("^MT-", gene),
                  rps_genes = grepl("^RPS", gene),
                  mrp_genes = grepl("^MRP", gene),
                  rpl_genes = grepl("^RPL", gene)) %>%
    dplyr::filter(mt_genes == FALSE &
                  rps_genes == FALSE &
                  mrp_genes == FALSE &
                  rpl_genes == FALSE &
                  t_value > 0) %>%
    arrange(state, -t_value) %>%
    group_by(state) %>%
    dplyr::slice(1:200) %>%
    dplyr::select(state, gene) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~ .x[[1]])) %>%
    deframe()
}

# estimate state list
param_df <- param_df %>% 
  dplyr::mutate(gene_set = map(ct_folder, get_de_genes)) %>%
  dplyr::mutate(gene_set = map2(ct, gene_set, function(x, y) {
    gene_sets <- y
    names(gene_sets) <- paste0(x,"_", names(gene_sets))
    return(gene_sets)
  }))

saveRDS(param_df$gene_set %>% unlist(recursive = F), 
        file = "./results/ct_data/state_genesets.rds")


