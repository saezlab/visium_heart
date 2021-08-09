# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of pseudobulk profiles between atlases

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)
library(edgeR)
library(ComplexHeatmap)

source("./analysis/utils/pseudobulk_utils.R")

cor_fun = circlize::colorRamp2(c(-1, 0, 1), 
                               c("darkred", "white", "darkblue"))

dist_fun = circlize::colorRamp2(c(0, 1), 
                               c("darkblue", "white"))

# Define a function of cell comparison
self_compare <- function(expression_mat, jsd = T) {
  
  if(jsd) {
    compare_profiles(expression_mat,
                     expression_mat)
    
  } else {
    
    compare_profiles(expression_mat,
                     expression_mat,
                     jsd = F)
    
  }
  
}
  
# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_wstates.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi") %>%
  dplyr::mutate(filt_mats = map(value, edgeR_filtering, min.count = 5)) %>%
  dplyr::mutate(norm_mats = map(filt_mats, cpm_norm)) %>%
  dplyr::mutate(self_jsd = map(filt_mats, self_compare),
                self_cor = map(norm_mats, self_compare, jsd =F))

pdf("./visium_results_manuscript/hca/MI_self_jsd.pdf")

walk(mi_atlas_mats$self_jsd, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/MI_self_cor.pdf")

walk(mi_atlas_mats$self_cor, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()

# HC atlas, we could expect different grouping vars:
hca_atlas <- readRDS("./visium_results_manuscript/integration/hca_pseudobulk.rds")
hca_atlas_mats <- map(hca_atlas, function(x) {
  
  atlas_mats <- x[names(x) != "annotations"]
  atlas_mats <- map(atlas_mats, assay)  %>%
    enframe("ann_level")

}) %>%
  enframe("origin") %>%
  unnest() %>%
  mutate(atlas = "hca") %>%
  dplyr::mutate(filt_mats = map(value, edgeR_filtering, min.count = 5)) %>%
  dplyr::mutate(norm_mats = map(filt_mats, cpm_norm)) %>%
  dplyr::mutate(self_jsd = map(filt_mats, self_compare),
                self_cor = map(norm_mats, self_compare, jsd =F))

# Self reports
pdf("./visium_results_manuscript/hca/hca_self_jsd.pdf", width = 12, height = 12)

walk(hca_atlas_mats$self_jsd, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/hca_self_cor.pdf", width = 12, height = 12)

walk(hca_atlas_mats$self_cor, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()

# Compare everything versus everything

hca_atlas_mats <- hca_atlas_mats %>%
  dplyr::mutate(dist_mats_cell_type = map(value, 
                                compare_profiles, 
                                matB = mi_atlas_mats$value[[1]]),
                cor_mats_cell_type = map(norm_mats, 
                                          compare_profiles, 
                                          matB = mi_atlas_mats$norm_mats[[1]],
                                         jsd = F),
                dist_mats_cell_state = map(value, 
                                          compare_profiles, 
                                          matB = mi_atlas_mats$value[[2]]),
                cor_mats_cell_state = map(norm_mats, 
                                         compare_profiles, 
                                         matB = mi_atlas_mats$norm_mats[[2]],
                                         jsd = F))



pdf("./visium_results_manuscript/hca/js_distances_hca_celltypes.pdf")

walk(hca_atlas_mats$dist_mats_cell_type, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/js_distances_hca_cellstates.pdf",
    height = 20, width = 15)

walk(hca_atlas_mats$dist_mats_cell_state, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()


pdf("./visium_results_manuscript/hca/cor_distances_hca_celltypes.pdf")

walk(hca_atlas_mats$cor_mats_cell_type, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/cor_distances_hca_cellstates.pdf",
    height = 20, width = 15)

walk(hca_atlas_mats$cor_mats_cell_state, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()


