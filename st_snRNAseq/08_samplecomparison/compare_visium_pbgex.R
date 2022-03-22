# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we analyse the pseudobulk profiles of visium data sets
#' We go full supervised from the get_patient_annotations_revisions grouping (We don't perform clustering)

library(tidyverse)
library(scater)
library(uwot)
library(cowplot)
library(ggrepel)
library(factoextra)
source("./analysis/utils/pseudobulk_utils.R")

# Patient annotations ------------------------------------------------------
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
pseudobulk_data <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["gex"]]
#atlas_meta <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["annotations"]]

pb_meta <- colData(pseudobulk_data) %>% 
  as.data.frame() %>% 
  left_join(sample_dict, 
            by = c("slide.meta.data...vars." = "sample_id"))

patient_cells <- pb_meta %>%
  group_by(patient_id) %>%
  summarize(ncells = sum(ncells))

# This summarizes the info into patients
pb_data <- sumCountsAcrossCells(assay(pseudobulk_data), 
                                DataFrame(pb_meta[, c("patient_id")]))


# Normalizing and filtering genes ------------------------------------------------------
pseudobulk_data_counts <- assay(pb_data)
colnames(pseudobulk_data_counts) <- pb_data$pb_meta...c..patient_id...
pseudobulk_data_counts <- edgeR_filtering(pseudobulk_data_counts, 
                                          min.count = 100,
                                          min.prop = 0.85,
                                          min.total.count = 100)

norm_pseudobulk_data_counts <- cpm_norm(pseudobulk_data_counts)

# UMAP ------------------------------------------------------
set.seed(241099)

gex_umap <- umap(t(norm_pseudobulk_data_counts), 
                 n_neighbors = 7, 
                 n_epochs = 1000,
                 metric = "cosine") %>%
  as.data.frame() %>%
  mutate(patient_id = colnames(norm_pseudobulk_data_counts))

gex_umap_dat <- gex_umap %>%
  left_join(sample_dict %>% dplyr::select(patient_id, patient, patient_group) %>% unique())

write.table(gex_umap_dat, 
            col.names = T, row.names = F, quote = F, 
            sep = ",", file = "./results/sample_comparison/visium_pb/visium_pb_umap.txt")

pdf("./results/sample_comparison/visium_pb/visium_pb_umap.pdf", height = 5, width = 6)

plot(gex_umap_dat %>%
ggplot(aes(x = V1, y = V2, 
           color = patient, 
           label = patient_id)) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  xlab("UMAP1") +
  ylab("UMAP2"))

plot(gex_umap_dat %>%
  ggplot(aes(x = V1, y = V2, 
             color = patient_group, 
             label = patient_id)) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  xlab("UMAP1") +
  ylab("UMAP2") )

dev.off()

# A better quantitative version is clustering

gex_hclust <- eclust(t(norm_pseudobulk_data_counts), "hclust", k = 3)

color_palette <- tibble(patient_id = gex_hclust$labels[gex_hclust$order]) %>%
  left_join(sample_dict[,c("patient_group", "patient_id")] %>% unique()) %>%
  left_join(tibble(patient_group = c("group_1", "group_2", "group_3"),
                   col = c("red", "darkgreen", "blue")))

pdf("./results/sample_comparison/visium_pb/visium_pb_hclust.pdf", height = 8, width = 5)

plot(fviz_dend(gex_hclust, 
               rect = TRUE, 
               label_cols = color_palette$col,
               k_colors = rep("black",3)))

dev.off()

pdf("./results/sample_comparison/visium_pb/visium_pb_umap_v2.pdf", height = 3.5, width = 4)
plt <- read_csv("./results/sample_comparison/visium_pb/visium_pb_umap.txt") %>%
  ggplot(aes(x = V1, y = V2, color = patient_group)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none") +
  ylab("UMAP1 (visium pseudobulk)") +
  xlab("UMAP2 (visium pseudobulk)")

plot(plt)
dev.off()
