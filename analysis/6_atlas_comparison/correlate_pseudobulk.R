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

# Dictionary between atlases
ct_dictionary <- read.table("./markers/hca_mi_ctannotations.txt", sep = "\t", header = T)

# Aesthetics definitions

cor_fun = circlize::colorRamp2(c(-1, 0, 1), 
                               c("darkred", "white", "darkblue"))

dist_fun = circlize::colorRamp2(c(0, 1), 
                               c("darkblue", "white"))


#' Compare expression profiles of the same matrix
#' @param expression_mat: a expression matrix with genes in rows and samples in colums
#' @param jsd: do you want to use jensen shannon divergences distances?
#' @return a distance matrix
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

#' Filter matrix by gene
#' @param expression_mat: a expression matrix with genes in rows and samples in colums
#' @param gene_list: do you want to use jensen shannon divergences distances?
#' @return a reduced matrix containing only genes of interest
filter_genes <- function(expression_mat, gene_list) {
  
  gene_list <- gene_list[gene_list %in% rownames(expression_mat)]
  
  if(length(gene_list) > 0) {
    
    return(expression_mat[gene_list,])
    
  } else{
    
    NULL
    
  }
}

# Defining a set of genes that's useful for comparison

cell_type_mrkrs <- (read.table("./markers/Kuppe_def.txt",
                                                 sep = "\t",stringsAsFactors = F,
                                                 header = T))[,c(1,2)]

cell_type_mrkrs[,1] <- toupper(cell_type_mrkrs[,1])
colnames(cell_type_mrkrs) <- c("gene","cell_type")

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_fordeconv.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi") %>%
  dplyr::mutate(filt_mats = map(value, edgeR_filtering, min.count = 5)) %>%
  dplyr::mutate(norm_mats = map(filt_mats, cpm_norm)) %>%
  dplyr::mutate(filt_mats = map(filt_mats, filter_genes, gene_list = cell_type_mrkrs$gene %>% unique())) %>%
  dplyr::mutate(norm_mats = map(norm_mats, filter_genes, gene_list = cell_type_mrkrs$gene %>% unique())) %>%
  dplyr::mutate(self_jsd = map(filt_mats, self_compare),
                self_cor = map(norm_mats, self_compare, jsd =F))

pdf("./visium_results_manuscript/hca/MI_self_reduced_jsd.pdf")

walk(mi_atlas_mats$self_jsd, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/MI_self__reduced_cor.pdf")

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
  dplyr::mutate(filt_mats = map(filt_mats, filter_genes, gene_list = cell_type_mrkrs$gene %>% unique())) %>%
  dplyr::mutate(norm_mats = map(norm_mats, filter_genes, gene_list = cell_type_mrkrs$gene %>% unique())) %>%
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

# Paper comparisons

cell_type_nuclei_cor <- hca_atlas_mats %>% dplyr::filter(origin == "nuclei",
                                                         ann_level == "cell_type") %>%
  pull(cor_mats_cell_type)

cell_state_nuclei_cor <- hca_atlas_mats %>% dplyr::filter(origin == "nuclei",
                                                          ann_level == "cell_state") %>%
  pull(cor_mats_cell_state)

# Clean the matrices
cell_type_nuclei_cor[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column("HCA") %>%
  pivot_longer(-HCA, names_to = "MI") %>%
  mutate(HCA = gsub("[A-Z]_","",HCA),
         MI = gsub("[A-Z]_","", MI)) %>%
  left_join(ct_dictionary, by = c("HCA" = "HCA_ann")) %>%
  na.omit() %>%
  mutate(HCA = MI_ann) %>%
  ggplot(aes(x = HCA, 
             y = MI,
             fill = value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.1),
        axis.text = element_text(size = 12)) +
  scale_fill_gradient(
    low = "black",
    high = "yellow")
  

pdf("./visium_results_manuscript/hca/js_distances_hca_celltypes_reduced.pdf")

walk(hca_atlas_mats$dist_mats_cell_type, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/js_distances_hca_cellstates_reduced.pdf",
    height = 20, width = 15)

walk(hca_atlas_mats$dist_mats_cell_state, function(x) {
  
  draw(Heatmap(x, col = dist_fun))
  
} )

dev.off()


pdf("./visium_results_manuscript/hca/cor_distances_hca_celltypes_reduced.pdf")

walk(hca_atlas_mats$cor_mats_cell_type, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()

pdf("./visium_results_manuscript/hca/cor_distances_hca_cellstates_reduced.pdf",
    height = 20, width = 15)

walk(hca_atlas_mats$cor_mats_cell_state, function(x) {
  
  draw(Heatmap(x, col = cor_fun))
  
} )

dev.off()


