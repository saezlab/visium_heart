# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we compare the compositions of complete atlases
#' MI:
#' RNA, ATAC, Visium
#' 
#' First we show that our atlas is consistent and that per-patient, the correlation is quite high
#' 
#' 

library(tidyverse)

# Read individual compositions

atac_props <- read.csv("./results/compositions/atac_compositions.txt", sep = "\t")
rna_props <- read.csv("./results/compositions/snrna_compositions.txt", sep = "\t")
spatial_props <- read.csv("./results/compositions/spatial_compositions.txt", sep = "\t")

apex_samples <- c("Visium_11_CK289",
                  "Visium_12_CK290",
                  "Visium_13_CK291",
                  "Visium_14_CK292",
                  "Visium_15_CK293",
                  "Visium_16_CK294",
                  "Visium_18_CK296",
                  "Visium_20_CK298")


# patient grouping
pat_anns <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  dplyr::mutate(apex = ifelse(sample_id %in% apex_samples, "LV_APEX", "LV")) %>%
  dplyr::select(patient_id, patient_group, apex) %>%
  unique()

# All proportions
mi_props <- left_join(spatial_props, rna_props, 
                      by = c("patient_id", "cell_type")) %>%
  left_join(atac_props,
            by = c("patient_id", "cell_type")) %>%
  mutate(sn_n_cells = ifelse(is.na(sn_n_cells), 0, sn_n_cells),
         sn_prop_cells = ifelse(is.na(sn_prop_cells), 0, sn_prop_cells),
         atac_n_cells = ifelse(is.na(atac_n_cells), 0, atac_n_cells),
         atac_prop_cells = ifelse(is.na(atac_prop_cells), 0, atac_prop_cells)) %>%
  left_join(pat_anns)

# Here we take the mean estimated proportion per patient from all views
multiview_props <- mi_props %>%
  pivot_longer(-c("patient_id", "cell_type", "patient_group","apex")) %>%
  dplyr::filter(grepl("prop", name)) %>%
  group_by(patient_id, cell_type, apex) %>%
  summarise(multiview_mean_prop = mean(value)) %>%
  left_join(pat_anns) %>%
  dplyr::filter(patient_group != "group_1")

# From here we can estimate how many APEX samples we have per group
pat_plot <- multiview_props %>%
  ungroup() %>%
  dplyr::select(patient_id, apex, patient_group) %>%
  unique() %>%
  group_by(patient_group, apex) %>%
  summarize(npats = length(patient_id)) %>%
  ggplot(aes(x = apex,
             y = npats)) +
  geom_bar(stat = "identity") +
  facet_grid(.~patient_group)

pdf("./results/compositions/apex_overview.pdf", height = 4, width = 4)

plot(pat_plot)

dev.off()

# Now, we will correlate the compositions from the LV and APEX within a group

apex_cors <- multiview_props %>%
  group_by(apex, patient_group, cell_type) %>%
  summarize(multiview_mean_prop = median(multiview_mean_prop)) %>%
  group_by(patient_group) %>%
  nest() %>%
  mutate(data = map(data, function(dat) {
    pivot_wider(dat, names_from = "apex", values_from = "multiview_mean_prop")
  })) %>%
  mutate(cor_res = map(data, function(dat) {
    
    cor.test(dat$LV, dat$LV_APEX, method = "spearman") %>%
      broom::tidy()
    
  }))


apex_cors %>% 
  select(cor_res) %>%
  unnest() %>%
  dplyr::select(patient_group, estimate)
  

pdf("./results/compositions/apex_correlation.pdf", height = 3, width = 7)

plot(apex_cors %>%
  select(data) %>%
  unnest() %>%
  ggplot(aes(x = log1p(LV), y = log1p(LV_APEX), color = cell_type)) +
  geom_point() +
  xlim(c(0,.4)) +
  ylim(0,0.4) +
  theme_classic() +
  facet_grid(. ~ patient_group))

dev.off()


  




