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
sample_dict <- readRDS("./markers/snrna_patient_anns_revisions.rds")

excluded_samples <- sample_dict %>%
  dplyr::filter(patient_group != "group_1") %>%
  pull(patient_id)

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]]

mi_atlas_meta <- mi_atlas[names(mi_atlas) == "annotations"][[1]] %>%
  dplyr::select(orig.ident, cell_type) %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id")) %>%
  dplyr::select(patient_id, cell_type) %>%
  group_by(patient_id, cell_type) %>%
  dplyr::summarize(ncells_ct = length(cell_type)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(ncells_sample = sum(ncells_ct)) %>%
  dplyr::mutate(propcells_ct = ncells_ct/ncells_sample) %>%
  dplyr::left_join(ct_dictionary, by = c("cell_type" = "MI_ann")) %>%
  dplyr::mutate(unified_ct = cell_type,
                sample = patient_id,
                atlas = "MI") %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::select(sample, unified_ct, propcells_ct, atlas)

# HC atlas, we could expect different grouping vars:
hca_atlas_meta <- readRDS("./ext_data/hca_pseudobulk.rds")[["nuclei"]][["annotations"]]  %>%
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

# First correlations between cell-types median compositions --------------------------

all_medians <- all_compositions %>% group_by(unified_ct, atlas) %>%
  summarize(median_prop = median(propcells_ct)) %>%
  pivot_wider(names_from = atlas, values_from = median_prop)

write_csv(all_medians, file = "./results/hca_comparison/hcacomp_med_prop_all.csv")

cor_all <- cor.test(log10(all_medians$HCA), log10(all_medians$MI), method = "spearman")

pdf("./results/hca_comparison/hcacomp_medprop_all.pdf", height = 3, width = 3)

ggplot(all_medians, aes(x = log10(HCA), y = log10(MI), label = unified_ct)) +
  ggrepel::geom_text_repel() +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ggtitle(paste0("Spearman Correlation, \n", 
                 round(cor_all$estimate, 3), 
                 ", p-value = ",
                 round(cor_all$p.value,3)))

dev.off()

# Then correlations between cell-types median compositions without disease pats
# This should look better --------------------------

healthy_medians <- all_compositions %>% 
  dplyr::filter(! sample %in% excluded_samples) %>%
  group_by(unified_ct, atlas) %>%
  summarize(median_prop = median(propcells_ct)) %>%
  pivot_wider(names_from = atlas, values_from = median_prop)

write_csv(healthy_medians,
            file = "./results/hca_comparison/hcacom_medprop_healthy.csv")

cor_healthy <- cor.test(log10(healthy_medians$HCA), log10(healthy_medians$MI), method = "spearman")

pdf("./results/hca_comparison/hcacom_medprop_healthy.pdf", height = 3, width = 3)

ggplot(healthy_medians, aes(x = log10(HCA), y = log10(MI), label = unified_ct)) +
  ggrepel::geom_text_repel() +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ggtitle(paste0("Spearman Correlation, \n", 
                 round(cor_healthy$estimate, 3), 
                 ", p-value = ",
                 round(cor_healthy$p.value,3)))

dev.off()

# Comparison of compositions ---------------------------------------------------------------

vln_plts <- ggplot(all_compositions, aes(y = propcells_ct, 
                             x = atlas,
                             fill = atlas)) + 
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = .05)) +
  facet_wrap(~unified_ct, ncol = 4,scales = "free_y") + theme_classic() +
  ylab("Proportions")


write.table(all_compositions, file = "./results/hca_comparison/all_compositions.txt", col.names = T, row.names = F, quote = F, )

pdf("./results/hca_comparison/compositions_box_all.pdf")

plot(vln_plts)

dev.off()

# Comparison of compositions (only healthy) --------------------------------------------------------

healthy_compositions <- all_compositions %>%
  dplyr::filter(! sample %in% excluded_samples)
  
vln_plts_h <- ggplot(healthy_compositions, aes(y = propcells_ct, 
                                         x = atlas,
                                         fill = atlas)) + 
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = .05)) +
  facet_wrap(~unified_ct, ncol = 4,scales = "free_y") + theme_classic() +
  ylab("Proportions")


write.table(healthy_compositions, file = "./results/hca_comparison/healthy_compositions.txt", col.names = T, row.names = F, quote = F, )

pdf("./results/hca_comparison/compositions_box_healthy.pdf")

plot(vln_plts_h)

dev.off()

# Proper statistical tests

# First wilcoxon tests between shared cell-types

wilcox_all <- all_compositions %>%
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

write.table(wilcox_all, file = "./results/hca_comparison/wilcox_all.txt", col.names = T, row.names = F, quote = F)


wilcox_healthy <- healthy_compositions %>%
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

write.table(wilcox_healthy, file = "./results/hca_comparison/wilcox_healthy.txt", col.names = T, row.names = F, quote = F)

