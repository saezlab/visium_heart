# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of pseudobulk profiles between atlases

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)

norm_mat <- function(x, scale_factor = 10000) {
  norm_mtx <- log1p(((x + 1)/colSums(x)) * scale_factor)
}

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_annotated.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi")
  
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
    
    dist_mat <- cor((norm_mat(cbind(matA, matB))))
    dist_mat <- dist_mat[colnames(matA),colnames(matB)]
    
    return(dist_mat)
  }
}

# Compare everything versus everything

hca_atlas_mats <- hca_atlas_mats %>%
  dplyr::mutate(dist_mats = map(value, 
                                compare_profiles, 
                                matB = mi_atlas_mats$value[[1]]))

pdf("./visium_results_manuscript/hca/js_distances_hca.pdf")

walk(hca_atlas_mats$dist_mats, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()

