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
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_wstates.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi")

# Annotations to check I didn't mess up
mi_anns <- mi_atlas["annotations"]$annotations
cts <- unique(mi_anns$cell_type)

map(set_names(cts), function(x) {
  
  filt_anns <- dplyr::filter(mi_anns,
                             cell_type == x)
  
  ncells <- nrow(filt_anns)
  
  nstates <- sum(grepl(x, filt_anns$deconv_col))
  
  ncells == nstates
})


# I will look at the correlation of my states
comp_mat <- 
self_corr <- JSD(t(mi_atlas_mats$value[[2]]),est.prob = "empirical") ** (1/2)
colnames(self_corr) = rownames(self_corr) <- colnames(mi_atlas_mats$value[[2]])
corrplot(self_corr, is.corr = F,method = "color")

mds_fit <- as.data.frame(cmdscale(as.dist(sample_divergences), 
                                  eig=TRUE, k=2)$points) %>%
  rownames_to_column("snRNA_state")



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
  dplyr::mutate(dist_mats_cell_type = map(value, 
                                compare_profiles, 
                                matB = mi_atlas_mats$value[[1]]),
                dist_mats_cell_state = map(value, 
                                          compare_profiles, 
                                          matB = mi_atlas_mats$value[[2]]))

pdf("./visium_results_manuscript/hca/js_distances_hca_celltypes.pdf")

walk(hca_atlas_mats$dist_mats_cell_type, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()

pdf("./visium_results_manuscript/hca/js_distances_hca_cellstates.pdf")

walk(hca_atlas_mats$dist_mats_cell_state, corrplot, 
     is.corr = FALSE, method = "color", 
     tl.cex = 0.6, tl.col = "black")

dev.off()
