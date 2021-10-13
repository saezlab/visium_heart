# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we compare the compositions of complete atlases
#' As supplementary, we also correlate patient specific compositions

library(tidyverse)


# ATAC DATA -----------------------------------------------------------------

atac_patient_info <- readRDS("./markers/atac_patient_anns_revisions.rds") 

atac_atlas_meta <- read_table2("./processed_atac/snATAC_annotated_meta.txt") %>%
  left_join(atac_patient_info, by = c("orig.ident" = "sample_id")) %>%
  dplyr::select(cell_type, patient_id)

atac_compositions <- atac_atlas_meta %>%
  group_by(patient_id, cell_type) %>%
  dplyr::summarize(ncells_ct = length(cell_type)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(ncells_sample = sum(ncells_ct)) %>%
  dplyr::mutate(propcells_ct = ncells_ct/ncells_sample) %>%
  dplyr::mutate(atlas = "ATAC") %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "Pericyte", "PC", cell_type)) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "neuronal", "Neuronal", cell_type)) 

# RNAseq DATA -----------------------------------------------------------------

# Dictionary between atlases
ct_dictionary <- read.table("./markers/hca_mi_ctannotations.txt", sep = "\t", header = T)
sample_dict <- readRDS("./markers/snrna_patient_anns_revisions.rds")

excluded_samples <- sample_dict %>%
  dplyr::filter(patient_group != "group_1") %>%
  pull(patient_id)

# Our atlas, we could expect different grouping vars:
mi_atlas <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]]

rna_compositions <- mi_atlas[names(mi_atlas) == "annotations"][[1]] %>%
  dplyr::select(orig.ident, cell_type) %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id")) %>%
  dplyr::select(patient_id, cell_type) %>%
  group_by(patient_id, cell_type) %>%
  dplyr::summarize(ncells_ct = length(cell_type)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(ncells_sample = sum(ncells_ct)) %>%
  dplyr::mutate(propcells_ct = ncells_ct/ncells_sample) %>%
  dplyr::mutate(atlas = "RNA") %>%
  ungroup()

# All data

all_props <- bind_rows(atac_compositions, rna_compositions)

median_props <- all_props %>%
  dplyr::group_by(cell_type, atlas) %>%
  summarize(median_prop = median(propcells_ct)) %>%
  pivot_wider(names_from = atlas, values_from = median_prop, values_fill = 0)

cor_all <- cor.test(log1p(median_props$RNA), log1p(median_props$ATAC), method = "spearman")


ggplot(median_props, aes(x = log1p(RNA), y = log1p(ATAC), label = cell_type, color = cell_type)) +
  ggrepel::geom_text_repel() +
  geom_point() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  ggtitle(paste0("Spearman Correlation, ", 
                 round(cor_all$estimate, 3), 
                 ", p-value = ",
                 round(cor_all$p.value,3)))


patient_props <- all_props %>%
  dplyr::select(-c("ncells_ct", "ncells_sample")) %>%
  group_by(patient_id) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    
    x %>%
      dplyr::group_by(cell_type, atlas) %>%
      summarize(median_prop = median(propcells_ct)) %>%
      pivot_wider(names_from = atlas, values_from = median_prop, values_fill = 0)
    
  })) %>%
  dplyr::filter(patient_id != "IZ_P11")


patient_props %>%
  mutate(cor_results = map(data, function(dat) {
    
    cor.test(log1p(dat$RNA), log1p(dat$ATAC), method = "spearman") %>%
      broom::tidy()
    
    
  })) %>%
  select(cor_results) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value))



