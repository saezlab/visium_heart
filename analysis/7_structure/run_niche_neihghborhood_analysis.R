# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We perform a neighborhood analysis inside thee
#' Giotto universe

library(Giotto)
library(Seurat)
library(tidyverse)

# Defining python path for Giotto
my_python_path <- "/Users/ricardoramirez/opt/miniconda3/envs/spat_giotto/bin/python"

instrs <- createGiottoInstructions(python_path = my_python_path,
                                  show_plot = TRUE,  
                                  return_plot = TRUE,
                                  save_plot = FALSE,
                                  plot_format = 'png',
                                  dpi = 200,
                                  height = 9,
                                  width = 9)

# These are the niches that define the structures in many slides
selected_niches <- readRDS("./markers/main_niches.rds")

# Define df of parameters
path <- "./visium_data/"

processed_path <- "./processed_visium/objects/"

sample_names <- list.files(path)

slide_files <- paste0(path,
                      sample_names,
                      "/outs")

seurat_files <- paste0(processed_path,
                       sample_names,
                       ".rds")

param_df <- tibble(sample_name = sample_names,
                   slide_file = slide_files,
                   seurat_file = seurat_files,
                   enrichment_out = paste0("./results/niche_mapping/giotto/result_tables/",
                                           sample_names, "_neighbor_res.rds"),
                   network_out = paste0("./results/niche_mapping/giotto/plots/",
                                        sample_names, "_neighbor_res.pdf"))

run_neighborhood_analysis <- function(sample_name, slide_file, seurat_file, enrichment_out, network_out) {
  
  print(sample_name)
  
  # Reading Giotto from raw data
  giotto_obj <- Giotto::createGiottoVisiumObject(h5_visium_path = paste0(slide_file, "/filtered_feature_bc_matrix.h5"),
                                                 h5_tissue_positions_path = paste0(slide_file, "/spatial/tissue_positions_list.csv"),
                                                 instructions = instrs)
  # Reading Seurat from processed data
  seurat_obj <- readRDS(seurat_file)
  
  # Getting the data of meaningful niches
  niche_meta <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, opt_clust_integrated) %>%
    group_by(opt_clust_integrated) %>%
    summarise(n_spots = length(opt_clust_integrated)) %>%
    mutate(prop_spots = n_spots/sum(n_spots)) %>%
    dplyr::filter(prop_spots > 0.01)
  
  slide_selected_niches <- selected_niches[selected_niches %in% niche_meta$opt_clust_integrated]
  
  slide_meta <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(opt_clust_integrated %in% slide_selected_niches) %>%
    column_to_rownames("cell_id")
  
  # Subsetting the object
  giotto_obj <- subsetGiotto(giotto_obj, 
                             cell_ids = rownames(slide_meta))
  
  giotto_obj@cell_metadata$niche <- factor(slide_meta[giotto_obj@cell_metadata$cell_ID, 
                                                      "opt_clust_integrated"], 
                                           levels = sort(selected_niches))
  # Generate spatial network
  giotto_obj <- createSpatialNetwork(gobject = giotto_obj, k = 50)
  
  cell_proximities <- cellProximityEnrichment(gobject = giotto_obj, 
                                             cluster_column = 'niche',
                                             spatial_network_name = 'Delaunay_network',
                                             adjust_method = 'fdr',
                                             number_of_simulations = 1000)
  #barplot
  net_plot <- cellProximityNetwork(gobject = giotto_obj,
                                   CPscore = cell_proximities,
                                   remove_self_edges = T, 
                                   only_show_enrichment_edges = T)
  
  pdf(file = network_out)
  
  print(net_plot)
  
  dev.off()
  
  saveRDS(cell_proximities$enrichm_res %>% as.data.frame(), file = enrichment_out)
}

# Main
pwalk(param_df, run_neighborhood_analysis)

