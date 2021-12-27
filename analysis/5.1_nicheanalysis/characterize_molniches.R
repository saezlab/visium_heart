# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run characterization of molecular niche definitions 

library(Seurat)
library(tidyverse)
source("./analysis/utils/niche_utils.R")

spot_anns <- readRDS("./processed_visium/integration/integrated_slides_umap.rds")[["meta_data"]] %>%
  as.data.frame() %>%
  rownames_to_column("raw_spot_id") %>%
  dplyr::mutate(spot_id = strsplit(raw_spot_id, "_") %>%
                  map_chr(., ~ .x[[1]]) %>%
                  paste0(orig.ident, "..", .))

umap_data <- readRDS("./processed_visium/integration/integrated_slides_umap.rds")[["reduction"]]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("raw_spot_id") %>%
  left_join(spot_anns, by = "raw_spot_id") %>%
  dplyr::select(spot_id, UMAP_1, UMAP_2)


pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

# Select resolutions of interest

int_res <- c("Spatial_snn_res.0.2", "Spatial_snn_res.0.4" , "Spatial_snn_res.0.5", "Spatial_snn_res.0.6", "Spatial_snn_res.0.8")

spot_anns <- spot_anns %>%
  left_join(pat_anns, by = c("orig.ident" = "sample_id")) %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", int_res))

ct_description <- readRDS("./results/ind_mats/cell_type_props.rds")
state_description <- readRDS("./results/ind_mats/cell_state_scorespos.rds")
  
#sample_spot

for(res in int_res) {
  
  print(res)
  
  # Create directory
  
  res_dir <- paste0("./results/niche_mapping/", res)
  
  system(paste0("mkdir ", res_dir))
  
  # Get the information of resolution of interest
  
  niche_info <- spot_anns %>%
    select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
    rename("mol_niche" = res) %>%
    dplyr::mutate(mol_niche = paste0("niche_", mol_niche))
  
  # 0. Plot UMAP
  
  umap_plt <- niche_info %>%
    left_join(umap_data) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = mol_niche)) +
    theme_classic() +
    ggrastr::geom_point_rast(size = 0.3) +
    guides(colour = guide_legend(override.aes = list(size=3)))

  pdf(file = paste0(res_dir, "/umap.pdf"), height = 5, width = 6)
  
  plot(umap_plt)
  
  write_csv(niche_info %>%
              left_join(umap_data), paste0(res_dir, "/umap.csv"))
  
  dev.off()
  
  # 1. Get the proportions of niche per sample id
  
  niche_props <- get_niche_props(niche_info)
  
  # 1.1 make dotplot
  
  pdf(file = paste0(res_dir, "/niche_props.pdf"), height = 7, width = 5)
  
  plot(plot_dots_niche(niche_props))
  
  write_csv(niche_props, paste0(res_dir, "/niche_props.csv"))
  
  dev.off()
  
  # 2. Filter patients that weren't described by niches
  
  filtered_niche_props <- filter_compositions(niche_props)
  
  # 3. Perform Kruskal-Wallis
  
  kw_niche_prop_test(filtered_niche_props = filtered_niche_props) %>%
    write_csv(paste0(res_dir, "/niche_props_kwtest.csv"))
  
  # 4. Plot boxplots
  
  pdf(file = paste0(res_dir, "/niche_patcomp.pdf"), height = 16, width = 12)
  
  plot(plot_box_niches(filtered_niche_props = filtered_niche_props))
  write_csv(filtered_niche_props, paste0(res_dir, "/niche_patcomp.csv"))
  
  dev.off()
  
  # 5. How much variance can be explained by the compositions of these niches?
  
  niche_prop_mat <- nichedf_tomatrix(filtered_niche_props)
  
  # ILR transform
  
  ILR_mat <- ILR_transform(niche_prop_mat = niche_prop_mat)
  
  # Explained variance
  explvar <- estimate_explvar(ILR_mat)
  
  explvar %>% write_csv(paste0(res_dir, "/explvar.csv"))
  
  explvar_val <- explvar %>%
    dplyr::filter(p.adj < 0.1) %>%
    pull(expl_var) %>%
    sum()
  
  # Clustering
  
  pdf(paste0(res_dir, "/niche_patclust.pdf"), height = 6, width = 5)
  
  plot(plot_clust(ILR_mat = ILR_mat,explvar_val = explvar_val))
  
  dev.off()
  
  # 6. What are the representative cell-types and cell-states
  
  run_typecharacterization(ct_description, niche_info, res_dir)
  
  run_statecharacterization(state_description = state_description,
                            niche_info, res_dir)
  
}







