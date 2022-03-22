# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate PROGENy scores with niches

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

# Main
spot_anns <- read_csv("./results/niche_mapping/mol_clust_class.csv")

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()

all_paths <- map(set_names(visium_df$visium_file, visium_df$sample), function(visium_file) { 
  print(visium_file)
  
  path_score <- readRDS(visium_file) %>%
    GetAssayData(., assay = "progeny") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spot_id") %>%
    pivot_longer(-spot_id)

  return(path_score)
})

all_paths <- enframe(all_paths) %>%
  dplyr::rename("orig.ident" = name) %>%
  unnest() %>%
  dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
  left_join(spot_anns) 

# For molecular

res <- "Spatial_snn_res.0.2"

niche_info <- all_paths %>%
  dplyr::filter(name != "TNFa") %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", res, "name", "value")) %>%
  dplyr::rename("mol_niche" = res) %>%
  dplyr::mutate(mol_niche = paste0("niche_", as.numeric(as.character(mol_niche)) + 1))

niche_info <- niche_info %>%
  group_by(mol_niche, name) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  group_by(name) %>%
  dplyr::mutate(std_value = (value - mean(value))/sd(value))

niche_info_mat <- niche_info %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = name, values_from = std_value) %>%
  column_to_rownames("mol_niche") %>%
  as.matrix()

act_plt <- ComplexHeatmap::Heatmap(t(niche_info_mat),name = "mean std act")

pdf(paste0("./results/niche_mapping/", res,"/progeny_characterization.pdf"), height = 4, width = 5)

ComplexHeatmap::draw(act_plt)

dev.off()

# For cell-type

res <- "composition_niche"

niche_info <- all_paths %>%
  dplyr::filter(name != "TNFa") %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", res, "name", "value")) %>%
  dplyr::rename("mol_niche" = res) %>%
  dplyr::mutate(mol_niche = paste0("niche_", mol_niche))

niche_info <- niche_info %>%
  group_by(mol_niche, name) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  group_by(name) %>%
  dplyr::mutate(std_value = (value - mean(value))/sd(value))

niche_info_mat <- niche_info %>%
  dplyr::select(-value) %>%
  pivot_wider(names_from = name, values_from = std_value) %>%
  column_to_rownames("mol_niche") %>%
  as.matrix()

act_plt <- ComplexHeatmap::Heatmap(t(niche_info_mat),name = "mean std act")

pdf(paste0("./results/niche_mapping/", res,"/progeny_characterization.pdf"), height = 4, width = 5)

ComplexHeatmap::draw(act_plt)

dev.off()

