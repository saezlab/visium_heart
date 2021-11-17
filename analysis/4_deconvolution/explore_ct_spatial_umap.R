# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Visualize marker genes in ILR umap of integrated visium data

library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)

sc_data <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")
sc_data <- scater::logNormCounts(sc_data)
umap_info <- readRDS("./results/niche_mapping/ct_niches/umap_compositional.rds")
cluster_info <- readRDS("./results/niche_mapping/ct_niches/niche_annotation_ct.rds")
integrated_compositions <- readRDS("./results/niche_mapping/ct_niches/integrated_compositions.rds")


# Add UMAP of ILR to the integrated object

meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  mutate(row_id = strsplit(row_id, "_") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(row_id = paste0(orig.ident, "..", row_id))

colnames(sc_data) <- meta_data$row_id

known_spots <- umap_info$row_id[umap_info$row_id %in% colnames(sc_data)]
sc_data <- sc_data[, known_spots]

umap_info <- umap_info %>%
  dplyr::filter(row_id %in% known_spots) %>%
  dplyr::select(V1, V2, row_id) %>%
  column_to_rownames("row_id") %>%
  as.matrix()

reducedDim(sc_data, type = "UMAP_ILR") <- umap_info[colnames(sc_data),] 

# Add niche labels

rownames(cluster_info) <- cluster_info$row_id
sc_data$niche <- cluster_info[colnames(sc_data), "ct_niche"]

# Add ct props

colData(sc_data) <- cbind(colData(sc_data), integrated_compositions[colnames(sc_data),])

pdf("./results/niche_mapping/ct_niches/ct_ILR_umap_exploration.pdf", height = 4, width = 4.5)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "niche")

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "CM", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Endo", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "vSMCs", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Lymphoid", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Myeloid", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "NPPB", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "TTN", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "MYH11", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "POSTN", point_size = 0.1)

dev.off()



