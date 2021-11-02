# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of pseudobulk profiles between atlases
#' We create long versions of cpm and raw counts and 
#' calculate correlations of cpms or distances with JSD
#' 
#' We limit the gene universe to the top 200 markers of the HCA

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)
library(edgeR)
library(ComplexHeatmap)

source("./analysis/utils/pseudobulk_utils.R")

# This is the ct naming convention dictionary -------------------------------------------
ct_dictionary <- read.table("./markers/hca_mi_ctannotations.txt", sep = "\t", header = T)

# Here we get the markers from the HCA --------------------------------------------------
nmarkers <- 200

hca_mrks <- read_csv("ext_data/hca_mrkrs.csv") %>%
  dplyr::filter(pvals_adj < 0.01) %>%
  arrange(group, -logfoldchanges) %>%
  dplyr::select(group, names) %>%
  group_by(group) %>%
  dplyr::slice(1:nmarkers) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

hca_mrks <- unlist(hca_mrks) %>% unique()
mi_mrks <- readRDS("./markers/pb_ct_marker_list.rds")
hca_mrks <- c(hca_mrks, unlist(mi_mrks)) %>% unique()

# We get our MI pb profile -------------------------------------------------------------------------
# This is pseudobulk profile per cell-type regardless of patient composition ------------------------
mi_atlas <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]]
mi_atlas_mat <- mi_atlas$gex
colnames_mat <- colData(mi_atlas$gex)[,1]
mi_atlas_mat <- assay(mi_atlas_mat)
colnames(mi_atlas_mat) <- colnames_mat

# We take genes with 0 everywhere ----------------------------------------------------------
mi_atlas_mat <- edgeR_filtering(mi_atlas_mat,
                                min.total.count = 100, 
                                min.count = 10, 
                                min.prop = 0.1)

# We calculate cpms -----------------------------------------------------------------------
mi_atlas_mat_cpm <- cpm_norm(mi_atlas_mat)
mi_atlas_mat_cpm_long <- mi_atlas_mat_cpm %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_type", values_to = "expr") %>%
  mutate(atlas = "MI")

# Let's repeat for HCA
hca_atlas <- readRDS("./ext_data/hca_pseudobulk.rds")[["nuclei"]]
hca_atlas <- assay(hca_atlas$cell_type)

# We take genes with 0 everywhere ----------------------------------------------------------
hca_atlas <- edgeR_filtering(hca_atlas,
                                min.total.count = 100, 
                                min.count = 10, 
                                min.prop = 0.1)

# We calculate cpms -----------------------------------------------------------------------
hca_atlas_cpm <- cpm_norm(hca_atlas)

hca_atlas_cpm_long <- hca_atlas_cpm %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell_type", values_to = "expr") %>%
  mutate(atlas = "HCA") %>%
  left_join(ct_dictionary, by = c("cell_type" = "HCA_ann")) %>%
  na.omit() %>%
  dplyr::select(-cell_type) %>%
  dplyr::rename("cell_type" = MI_ann)

# Merge all data

all_gex <- left_join(mi_atlas_mat_cpm_long, hca_atlas_cpm_long, by = c("cell_type", "gene")) %>%
  na.omit()

correlations_all <- all_gex %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(cor.res = map(data, function(d) {
    
    cor.test(d$expr.x, d$expr.y, method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  select(cor.res) %>%
  unnest()

correlations_filtered <- all_gex %>%
  dplyr::filter(gene %in% hca_mrks) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(cor.res = map(data, function(d) {
    
    cor.test(d$expr.x, d$expr.y, method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  select(cor.res) %>%
  unnest()

genes_all <- all_gex %>%
  pull(gene) %>%
  unique()

genes_filt <- all_gex %>%
  dplyr::filter(gene %in% hca_mrks) %>%
  pull(gene) %>%
  unique()

# Plots -------------------------------------------------------------------------------------
cells_filt <- ct_dictionary %>%
  na.omit()

mi_atlas_mat_cpm <- mi_atlas_mat_cpm[, cells_filt$MI_ann]
colnames(mi_atlas_mat_cpm) <- paste0("MI_", colnames(mi_atlas_mat_cpm) )

hca_atlas_cpm <- hca_atlas_cpm[, cells_filt$HCA_ann]
colnames(hca_atlas_cpm) <- cells_filt$MI_ann
colnames(hca_atlas_cpm) <- paste0("HCA_", colnames(hca_atlas_cpm))

cor_all <- cor(hca_atlas_cpm[genes_all,], mi_atlas_mat_cpm[genes_all,], method = "spearman")

cor_filt <- cor(hca_atlas_cpm[genes_filt,], mi_atlas_mat_cpm[genes_filt,], method = "spearman")

cor_all_tidy <- cor_all %>%
  as.data.frame() %>%
  rownames_to_column("atlas_a") %>%
  pivot_longer(-atlas_a, names_to = "atlas_b", values_to = "spearman_cor") %>%
  mutate(atlas_a = factor(atlas_a, 
                          levels = paste0("HCA_",cells_filt$MI_ann)),
         atlas_b =factor(atlas_b, 
                         levels = paste0("MI_",cells_filt$MI_ann)))

cor_filt_tidy <- cor_filt %>%
  as.data.frame() %>%
  rownames_to_column("atlas_a") %>%
  pivot_longer(-atlas_a, names_to = "atlas_b", values_to = "spearman_cor") %>%
  mutate(atlas_a = factor(atlas_a, 
                          levels = paste0("HCA_",cells_filt$MI_ann)),
         atlas_b =factor(atlas_b, 
                         levels = paste0("MI_",cells_filt$MI_ann)))

pdf("./results/hca_comparison/pseudobulk_exprcomp_.pdf", height = 4, width = 5)

plot(cor_filt_tidy %>%
  ggplot(aes(x = atlas_a, y = atlas_b, fill = spearman_cor)) +
  geom_tile() +
  coord_equal() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    ggtitle("Only marker genes")) #+
  #scale_fill_gradient2(midpoint = 0.5, limits = c(0,1))

plot(cor_all_tidy %>%
       ggplot(aes(x = atlas_a, y = atlas_b, fill = spearman_cor)) +
       geom_tile() +
       coord_equal() +
       theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust =0.5)) +
       ggtitle("All expression")) #+
  #scale_fill_gradient2(midpoint = 0.5, limits = c(0,1))

dev.off()
