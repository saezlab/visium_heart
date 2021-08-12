# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of abundances between atlases
#' We used HCA LV nuclei data for the sake of clarity in the text

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)
library(edgeR)
library(ComplexHeatmap)

source("./analysis/utils/pseudobulk_utils.R")

# Dictionary between atlases
ct_dictionary <- read.table("./markers/hca_mi_ctannotations.txt", sep = "\t", header = T)

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_fordeconv.rds")[[1]]
mi_atlas_meta <- mi_atlas[names(mi_atlas) == "annotations"][[1]] %>%
  dplyr::select(orig.ident, cell_type) %>%
  group_by(orig.ident, cell_type) %>%
  dplyr::summarize(ncells_ct = length(cell_type)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(ncells_sample = sum(ncells_ct)) %>%
  dplyr::mutate(propcells_ct = ncells_ct/ncells_sample) %>%
  dplyr::left_join(ct_dictionary, by = c("cell_type" = "MI_ann")) %>%
  dplyr::mutate(unified_ct = cell_type,
                sample = orig.ident,
                atlas = "MI") %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::select(sample, unified_ct, propcells_ct, atlas)


# HC atlas, we could expect different grouping vars:
hca_atlas_meta <- readRDS("./visium_results_manuscript/integration/hca_pseudobulk.rds")[["nuclei"]][["annotations"]]  %>%
  dplyr::select(sample, cell_type) %>%
  group_by(sample, cell_type) %>%
  dplyr::summarize(ncells_ct = length(cell_type)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(ncells_sample = sum(ncells_ct)) %>%
  dplyr::mutate(propcells_ct = ncells_ct/ncells_sample) %>%
  dplyr::left_join(ct_dictionary, by = c("cell_type" = "HCA_ann")) %>%
  dplyr::mutate(unified_ct = MI_ann,
                atlas = "HCA") %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::select(sample, unified_ct, propcells_ct, atlas)

# All data together ---------------------------------------------------------------

all_compositions <- bind_rows(mi_atlas_meta, hca_atlas_meta)

# Comparison of compositions ---------------------------------------------------------------

vln_plts <- ggplot(all_compositions, aes(y = propcells_ct, 
                             x = atlas,
                             fill = atlas)) + 
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = .05)) +
  facet_wrap(~unified_ct, ncol = 4,scales = "free_y") + theme_classic() +
  ylab("Proportions")

# Proper statistical tests

# First wilcoxon tests between shared cell-types

all_compositions %>%
  group_by(unified_ct) %>%
  nest() %>%
  mutate(wilcox_res = map(data, function(x) {
    wilcox.test(propcells_ct ~ atlas, 
                data = x, 
                alternative = "two.sided", 
                paired = F) %>% broom::tidy()
    
  })) %>%
  select(unified_ct, wilcox_res) %>%
  unnest()

# Second, stability of compositions

prop_variance <- all_compositions %>%
  group_by(unified_ct, atlas) %>%
  summarise(prop_sd = var(propcells_ct))

wilcox.test(prop_sd ~ atlas, 
            data = prop_variance, 
            alternative = "two.sided", 
            paired = F) %>% 
  broom::tidy()

variance_plot <- ggplot(prop_variance, aes(x = atlas, 
                                           y = prop_sd,
                                           fill = atlas)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = .05)) +
  ylab("std. dev. of cell-type proportions\
       between samples")  +
  theme(axis.text = element_text(size = "12"))


