# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we analyse the pseudobulk profiles of visium data sets
#' We go full supervised from the get_patient_annotations_revisions grouping (We don't perform clustering)

library(scater)
library(tidyverse)
library(philentropy)
library(uwot)
library(cowplot)
library(compositions)
library(ggrepel)
source("./analysis/utils/pseudobulk_utils.R")

# Patient annotations ------------------------------------------------------
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
pseudobulk_data <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["gex"]]
atlas_meta <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["annotations"]]

# Normalizing and filtering genes ------------------------------------------------------
pseudobulk_data_counts <- assay(pseudobulk_data)
colnames(pseudobulk_data_counts) <- pseudobulk_data$slide.meta.data...vars.
pseudobulk_data_counts <- edgeR_filtering(pseudobulk_data_counts, 
                                          min.count = 100,
                                          min.prop = 0.85,
                                          min.total.count = 100)

norm_pseudobulk_data_counts <- cpm_norm(pseudobulk_data_counts)

# UMAP ------------------------------------------------------
gex_umap <- umap(t(norm_pseudobulk_data_counts), 
                 n_neighbors = 5, 
                 n_epochs = 1000,
                 metric = "cosine") %>%
  as.data.frame() %>%
  mutate(orig.ident = colnames(norm_pseudobulk_data_counts))

gex_umap_dat <- gex_umap %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id"))

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
