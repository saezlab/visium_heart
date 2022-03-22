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
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "niche")

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "niche", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "niche", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "ident")
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "ident")

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "ident", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "ident", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "niche", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "orig.ident",  point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "CM", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "CM", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Endo", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Endo", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Fib", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Fib", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "vSMCs", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "vSMCs", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Lymphoid", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Lymphoid", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Myeloid", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Myeloid", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Adipo", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Adipo", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "PC", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "PC", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Mast", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Mast", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "prolif", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "prolif", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "Neuronal", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_HARMONY", colour_by = "Neuronal", point_size = 0.1)

plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "NPPB", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "TTN", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "MYH11", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "SPP1", point_size = 0.1)
plotReducedDim(sc_data, dimred = "UMAP_ILR", colour_by = "DCN", point_size = 0.1)

dev.off()

# Get genes that should be tested

# Add label to state
sc_data[["opt_state"]]  <- sc_data[["niche"]]

sc_meta <- colData(sc_data) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")
  

# Calculate the number of cells per state
# The minimum will represent the filtering of lowly expressed genes

print("Excluding lowly expressed genes")

n_cell_df <- colData(sc_data) %>% 
  as.data.frame() %>%
  group_by_at("opt_state") %>%
  summarise(n_spots = n())

min_spots <- n_cell_df %>%
  pull(n_spots) %>%
  min()

# First let's filter all genes that are extremely lowly expressed
expressed_genes <- counts(sc_data) > 0
gene_ix <- rowSums(expressed_genes) > 0
#gene_ix <- rowSums(expressed_genes) > min_spots
sc_data <- sc_data[gene_ix, ]
expressed_genes <- expressed_genes[gene_ix, ]

print("Keeping genes expressed in at least perc_thrsh in a state")

perc_thrsh = 0.5

# Then for each niche we will test specifically if the gene can be considered expressed
# Get cell_ids that belong to a class
niche_info <- colData(sc_data) %>% 
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  dplyr::select_at(c("cell_id", "opt_state")) %>%
  group_by_at("opt_state") %>%
  nest() %>%
  dplyr::rename("cell_ids" = data) %>%
  dplyr::mutate(cell_ids = map(cell_ids, ~ .x[[1]])) %>%
  left_join(n_cell_df) %>%
  dplyr::mutate(min_spots = (n_spots * perc_thrsh) %>% floor()) %>%
  dplyr::mutate(selected_genes = map2(cell_ids, min_spots, function(cids, mspots) {
    gene_ix <- rowSums(expressed_genes[, cids]) > mspots
    return(names(gene_ix[gene_ix]))
  }))

saveRDS(niche_info, file = "./results/niche_mapping/expressed_genes_niche.rds")


saveHDF5SummarizedExperiment(sc_data, dir = "./processed_visium/integration/integrated_slides_sce_ann/", 
                             prefix = "", replace = FALSE,
                             chunkdim = NULL, level = NULL, as.sparse = NA,
                             verbose = NA)

