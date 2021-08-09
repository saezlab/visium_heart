# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Compare previous annotations with current annotations

library(tidyverse)

previous_anns <- "./visium_results_manuscript/integration/scRNA.integrated_meta.rds"
current_anns <- "./visium_results_manuscript/integration/integrated_wstates_meta.rds"
umap_embedding <- "./visium_results_manuscript/integration/integrated_data_wstates_umap.rds"

previous_anns <- readRDS(previous_anns) %>%
  rownames_to_column("cell_id") %>%
  dplyr::mutate(cell_id = map_chr(strsplit(split = "_",x = cell_id), ~.x[[2]])) %>%
  dplyr::select(cell_id, orig.ident, top_annotation)

current_anns <- readRDS(current_anns) %>%
  rownames_to_column("raw_cell_id") %>%
  dplyr::mutate(cell_id = map_chr(strsplit(split = "-",x = raw_cell_id), ~.x[[1]])) %>%
  dplyr::select(raw_cell_id, cell_id, orig.ident, cell_type)
  
umap_embedding <- readRDS(umap_embedding) %>%
  as.data.frame() %>%
  rownames_to_column("raw_cell_id") %>%
  left_join(current_anns, by = "raw_cell_id") %>%
  left_join(previous_anns, by =c("cell_id", "orig.ident"))

pdf("./visium_results_manuscript/integration/comparison_past_anns.pdf", height = 7, width = 9)

print(ggplot(umap_embedding, aes(x = UMAP_1, 
                           y = UMAP_2, 
                           color = top_annotation)) +
  geom_point(size = 0.4, alpha = 0.7) + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=3))))

print(ggplot(umap_embedding, aes(x = UMAP_1, 
                           y = UMAP_2, 
                           color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.7) + 
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=3))))


dev.off()
