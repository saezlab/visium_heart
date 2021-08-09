# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of pseudobulk profiles between atlases (HCA vs manuscript version)

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)

norm_mat <- function(x, scale_factor = 10000) {
  norm_mtx <- log1p(((x + 1)/colSums(x)) * scale_factor)
}

# Function to get merged matrices from original atlas:
getmatrix <- function(mi_atlas_alt) {
  mi_atlas_alt <- mi_atlas_alt %>%
    map(function(x) { 
      as.data.frame(x) %>%
        rownames_to_column("gene") %>%
        pivot_longer(-gene)
    }) %>%
    enframe() %>% 
    dplyr::select(-name) %>%
    unnest() %>%
    group_by(name, gene) %>%
    summarise(raw_counts = sum(value)) %>%
    ungroup() %>%
    pivot_wider(values_from = raw_counts,
                names_from = gene,
                values_fill = 0)
  
  mi_atlas_cts <- mi_atlas_alt[["name"]]
  mi_atlas_alt <- as.matrix(mi_atlas_alt[,-1])
  rownames(mi_atlas_alt) <- mi_atlas_cts
  
  return(mi_atlas_alt)
}

# Implementation of distances
compare_profiles <- function(matA, matB, jsd = TRUE) {
  
  gene_ids <- intersect(rownames(matA), 
                        rownames(matB))
  
  colnames(matA) <- paste0("A_", colnames(matA))
  
  colnames(matB) <- paste0("B_", colnames(matB))
  
  matA <- matA[gene_ids, ]
  matB <- matB[gene_ids, ]
  
  if(jsd){
    
    dist_mat <- philentropy::JSD(t(cbind(matA, matB)), est.prob = "empirical")
    rownames(dist_mat) = colnames(dist_mat) <- c(colnames(matA),colnames(matB))
    dist_mat <- dist_mat ** (1/2)
    dist_mat <- dist_mat[colnames(matA),colnames(matB)]
    
    return(dist_mat)
    
  } else {
    
    dist_mat <- cor((norm_mat(cbind(matA, matB))),method = "spearman")
    dist_mat <- dist_mat[colnames(matA),colnames(matB)]
    
    return(dist_mat)
  }
}

# Our atlas previous version
mi_atlas <- readRDS("./visium_results_manuscript/pseudobulk/mi_pseudobulk.rds")
names(mi_atlas) <- strsplit(names(mi_atlas), "/") %>%
  map_chr(~ .x %>% last()) %>%
  strsplit("[.]") %>%
  map_chr(~ .x %>% first())

# Cell states
mi_atlas_alt <- map(mi_atlas, ~ assay(.x[["cell_type"]]))
mi_atlas_alt <- getmatrix(mi_atlas_alt)
mi_atlas_alt <- t(mi_atlas_alt)
mi_atlas_alt <- mi_atlas_alt[! rowSums(mi_atlas_alt) <= 100, ]

# Cell types
mi_atlas_mjr <- map(mi_atlas, ~ assay(.x[["major_cell_type"]]))
mi_atlas_mjr <- getmatrix(mi_atlas_mjr)
mi_atlas_mjr <- t(mi_atlas_mjr)
mi_atlas_mjr <- mi_atlas_mjr[! rowSums(mi_atlas_mjr) <= 100, ]

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_wstates.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi") %>%
  dplyr::mutate(dist_mats_cell_type = map(value, 
                                          compare_profiles, 
                                          matB = mi_atlas_alt),
                dist_mats_mjr = map(value, 
                                    compare_profiles, 
                                    matB = mi_atlas_mjr))

# HC atlas, we could expect different grouping vars:
hca_atlas <- readRDS("./visium_results_manuscript/integration/hca_pseudobulk.rds")
hca_atlas_mats <- map(hca_atlas, function(x) {
  
  atlas_mats <- x[names(x) != "annotations"]
  atlas_mats <- map(atlas_mats, assay)  %>%
    enframe("ann_level")
  
}) %>%
  enframe("origin") %>%
  unnest() %>%
  mutate(atlas = "hca")

hca_atlas_mats <- hca_atlas_mats %>%
  dplyr::mutate(dist_mats_cell_type = map(value, 
                                          compare_profiles, 
                                          matB = mi_atlas_alt),
                dist_mats_mjr = map(value, 
                                    compare_profiles, 
                                    matB = mi_atlas_mjr))

pdf("./visium_results_manuscript/hca/js_distances_hca_cellstates_old.pdf")

walk(hca_atlas_mats$dist_mats_cell_type, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()


pdf("./visium_results_manuscript/hca/js_distances_hca_celltypes_old.pdf")

walk(hca_atlas_mats$dist_mats_mjr, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()


pdf("./visium_results_manuscript/hca/js_distances_newMI_celltypes_old.pdf")

walk(mi_atlas_mats$dist_mats_mjr, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()


pdf("./visium_results_manuscript/hca/js_distances_newMI_cellstates_old.pdf")

walk(mi_atlas_mats$dist_mats_cell_type, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()


corrplot(compare_profiles(mi_atlas_mjr,mi_atlas_mjr,jsd = TRUE))

