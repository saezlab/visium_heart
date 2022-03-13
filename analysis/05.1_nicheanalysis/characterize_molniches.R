# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run characterization of molecular niche definitions 

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
source("./analysis/utils/niche_utils.R")

ct_niches <- readRDS("./results/niche_mapping/ct_niches/niche_annotation_ct.rds") %>%
  dplyr::select(row_id, niche) %>%
  dplyr::rename("spot_id" = row_id,
                "composition_niche" = niche)

spot_anns <- readRDS("./processed_visium/integration/integrated_slides_umap.rds")[["meta_data"]] %>%
  as.data.frame() %>%
  rownames_to_column("raw_spot_id") %>%
  dplyr::mutate(spot_id = strsplit(raw_spot_id, "_") %>%
                  map_chr(., ~ .x[[1]]) %>%
                  paste0(orig.ident, "..", .)) %>%
  left_join(ct_niches)

umap_data <- readRDS("./processed_visium/integration/integrated_slides_umap.rds")[["reduction"]]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("raw_spot_id") %>%
  left_join(spot_anns, by = "raw_spot_id") %>%
  dplyr::select(spot_id, UMAP_1, UMAP_2)

umap_info_ilr <- readRDS("./results/niche_mapping/ct_niches/umap_compositional.rds") %>%
  dplyr::select(V1, V2, row_id) %>%
  dplyr::rename("UMAP_1_ILR" = V1,
                "UMAP_2_ILR" = V2,
                "spot_id" = row_id)

pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

# Select resolutions of interest

int_res <- c(#"composition_niche",
             "Spatial_snn_res.0.2")#, "Spatial_snn_res.0.4" , 
             #"Spatial_snn_res.0.5","Spatial_snn_res.0.6", 
             #"Spatial_snn_res.0.8")

spot_anns <- spot_anns %>%
  left_join(pat_anns, by = c("orig.ident" = "sample_id")) %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", int_res))

ct_description <- readRDS("./results/ind_mats/cell_type_props.rds")
state_description <- readRDS("./results/ind_mats/cell_state_scorespos.rds")

#spot_anns %>%
#  mutate(raw_spot_id = strsplit(spot_id, "[.][.]") %>%
#           map_chr(., ~.x[[2]])) %>%
#  write_csv("./results/niche_mapping/mol_clust_class.csv")
  
#sample_spot

for(res in int_res) {
  
  print(res)
  
  # Create directory
  
  res_dir <- paste0("./results/niche_mapping/", res)
  
  system(paste0("mkdir ", res_dir))
  
  system(paste0("mkdir ", res_dir, "/niche_props_plts_patient"))
  system(paste0("mkdir ", res_dir, "/niche_props_plts_myogenic"))
  
  # Get the information of resolution of interest
  
  if(res == "composition_niche") {
    niche_info <- spot_anns %>%
      select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
      rename("mol_niche" = res) %>%
      dplyr::mutate(mol_niche = paste0("niche_", mol_niche))
    
  } else {
    niche_info <- spot_anns %>%
      select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
      rename("mol_niche" = res) %>%
      dplyr::mutate(mol_niche = paste0("niche_", as.numeric(as.character(mol_niche)) + 1))
  }
  
  
  
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
  
  # 0.1 Plot UMAP ILR
  
  umap_plt_ilr <- niche_info %>%
    left_join(umap_info_ilr) %>%
    ggplot(aes(x = UMAP_1_ILR, y = UMAP_2_ILR, color = mol_niche)) +
    theme_classic() +
    ggrastr::geom_point_rast(size = 0.3) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  pdf(file = paste0(res_dir, "/umap_ilr.pdf"), height = 5, width = 6)
  
  plot(umap_plt_ilr)
  
  write_csv(niche_info %>%
              left_join(umap_data), paste0(res_dir, "/umap_ilr.csv"))
  
  dev.off()
  
  # 1. Get the proportions of niche per sample id
  
  niche_props <- get_niche_props(niche_info)
  
  # 1.1 make dotplot
  
  pdf(file = paste0(res_dir, "/niche_props.pdf"), height = 7, width = 5)
  
  plot(plot_dots_niche(niche_props))
  
  write_csv(niche_props, paste0(res_dir, "/niche_props.csv"))
  
  dev.off()
  
  # 2.1 Filter niches that are patient specific
  
  filtered_niche_props <- filter_compositions(niche_props,by = "niche")
  
  # 2.2 Filter patients that weren't described by niches
  
  filtered_niche_props <- filter_compositions(filtered_niche_props, by = "patient")
  
  # 3. Perform Kruskal-Wallis
  
  kw_niche_prop_test(filtered_niche_props = filtered_niche_props) %>%
    write_csv(paste0(res_dir, "/niche_props_kwtest.csv"))
  
  kw_niche_prop_test_area(filtered_niche_props = filtered_niche_props) %>%
    write_csv(paste0(res_dir, "/niche_props_kwtest_area.csv"))
  
  kw_niche_prop_test_area(filtered_niche_props = filtered_niche_props, 
                          myogenic_enriched = T) %>%
    write_csv(paste0(res_dir, "/niche_props_kwtest_area_myogenic.csv"))
  
  # 3.1 Perform pairwise comparisons
  
  # Patient groups
  pw_niche_prop_comps <- pw_niche_prop_test(filtered_niche_props = filtered_niche_props,
                                            area = F,
                                            myogenic_enriched = F)
  
  pw_niche_prop_comps %>% write_csv(paste0(res_dir, "/niche_props_pwtest.csv"))
  
  # Myogenic enriched
  pw_niche_prop_comps_myog <- pw_niche_prop_test(filtered_niche_props = filtered_niche_props,
                     area = T
                     ,myogenic_enriched = T) 
  
  pw_niche_prop_comps_myog %>% write_csv(paste0(res_dir, "/niche_props_pwtest_area_myogenic.csv"))
  
  # 4. Plot boxplots
  
  # As individual plots
  pw_bplots_props <- plot_box_niches_pw(filtered_niche_props, pw_niche_prop_comps) %>%
    dplyr::mutate(pdf_out = paste0(res_dir, "/niche_props_plts_patient/", mol_niche, ".pdf"))
  
  walk2(pw_bplots_props$bxplt, pw_bplots_props$pdf_out, function(x,y) {
    
    pdf(y, height = 4, width = 2)
    
    plot(x)
    
    dev.off()
    
  })
  
  # As individual plots - myogenic 
  pw_bplots_props_myogenic <- plot_box_niches_pw(filtered_niche_props, pw_niche_prop_comps = pw_niche_prop_comps_myog, area = T) %>%
    dplyr::mutate(pdf_out = paste0(res_dir, "/niche_props_plts_myogenic/", mol_niche, ".pdf"))
  
  walk2(pw_bplots_props_myogenic$bxplt, pw_bplots_props_myogenic$pdf_out, function(x,y) {
    
    pdf(y, height = 4, width = 2)
    
    plot(x)
    
    dev.off()
    
  })
  
  
  pdf(file = paste0(res_dir, "/niche_patcomp.pdf"), height = 9, width = 12)
  
  plot(plot_box_niches(filtered_niche_props = filtered_niche_props))
  write_csv(filtered_niche_props, paste0(res_dir, "/niche_patcomp.csv"))
  
  dev.off()
  
  pdf(file = paste0(res_dir, "/niche_patcomp_area.pdf"), height = 16, width = 12)
  
  plot(plot_box_niches_area(filtered_niche_props = filtered_niche_props))
  write_csv(filtered_niche_props, paste0(res_dir, "/niche_patcomp_area.csv"))
  
  dev.off()
  
  pdf(file = paste0(res_dir, "/niche_patcomp_area_myogenic.pdf"), height = 9, width = 12)
  
  plot(plot_box_niches_area(filtered_niche_props = filtered_niche_props, 
                            myogenic_enriched = T))
  write_csv(filtered_niche_props, paste0(res_dir, "/niche_patcomp_area_myogenic.csv"))
  
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
  
  # Explained variance of individual regions
  
  explvar_area <- estimate_explvar_area(ILR_mat = ILR_mat, myogenic_enriched = F)
  
  explvar_area %>% write_csv(paste0(res_dir, "/explvar_area.csv"))
  
  explvar_area_val <- explvar_area %>%
    dplyr::filter(p.adj < 0.1) %>%
    pull(expl_var) %>%
    sum()
  
  explvar_area_me <- estimate_explvar_area(ILR_mat = ILR_mat, myogenic_enriched = T)
  
  explvar_area_me %>% write_csv(paste0(res_dir, "/explvar_area_me.csv"))
  
  explvar_area_me_val <- explvar_area_me %>%
    dplyr::filter(p.adj < 0.1) %>%
    pull(expl_var) %>%
    sum()
  
  # Explained variance of early time points
  
  explvar_early <- estimate_explvar_time(ILR_mat,early_only = T)
  
  explvar_early %>% write_csv(paste0(res_dir, "/explvar_time.csv"))
  
  explvar_val_early <- explvar_early %>%
    dplyr::filter(p.adj < 0.1) %>%
    pull(expl_var) %>%
    sum()
  
  # PCA of compositions of niches
  
  plot_PCA(ILR_mat = ILR_mat,
           explvar_val = explvar_val,
           early_only = FALSE,
           res_dir = res_dir)
  
  # PCA of composition of niches (Early time points)
  
  plot_PCA(ILR_mat = ILR_mat,
           explvar_val = explvar_val_early,
           early_only = T,
           res_dir = res_dir)
  
  plot_PCA_area(ILR_mat = ILR_mat,
           explvar_val = explvar_area_me_val,
           myogenic_enriched = T,
           res_dir = res_dir)
  
  plot_PCA_area(ILR_mat = ILR_mat,
                explvar_val = explvar_area_val,
                myogenic_enriched = F,
                res_dir = res_dir)
  
  
  # Clustering
  
  pdf(paste0(res_dir, "/niche_patclust.pdf"), height = 6, width = 5)
  
  plot(plot_clust(ILR_mat = ILR_mat,explvar_val = explvar_val))
  
  dev.off()
  
  
  #UMAP
  
  # 6. What are the representative cell-types and cell-states
  
  useful_niches <- unique(filtered_niche_props$mol_niche)
  niche_info_filt <- niche_info %>%
    dplyr::filter(mol_niche %in% useful_niches)
  
  
  run_characterization(ct_description = ct_description,
                       niche_info = niche_info,
                       res_dir = res_dir,
                       alias = "cell_type_enrich_wilcox",
                       height = 4,
                       width = 5)
  
  run_characterization(ct_description = state_description,
                       niche_info = niche_info,
                       res_dir = res_dir,
                       alias = "cell_state_enrich_wilcox",
                       height = 8,
                       width = 8)
  
  #run_typecharacterization(ct_description, niche_info, res_dir)
  
  #run_statecharacterization(state_description = state_description,
  #                          niche_info, res_dir)
  
}

# Generate a final annotation containing cell-states, cell-types and UMAP data

all_annotations <- spot_anns %>%
  mutate(raw_spot_id = strsplit(spot_id, "[.][.]") %>%
           map_chr(., ~.x[[2]])) %>%
  left_join(umap_data) %>%
  left_join(umap_info_ilr) %>%
  write_csv("./results/niche_mapping/niche_annotation_molct.csv")

# Add cell_type compositions (props)

integrated_compositions <- readRDS("./results/niche_mapping/ct_niches/integrated_compositions.rds") %>%
  as.data.frame() %>%
  rownames_to_column("spot_id")






