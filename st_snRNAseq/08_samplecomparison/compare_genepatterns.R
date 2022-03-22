# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare spark results between patient groups

library(tidyverse)

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")

pat_order <- sample_dict %>%
  dplyr::arrange(patient_group, patient_id) %>%
  pull(patient_id) %>%
  unique()

run_wilcox_all <- function(prop_data) {
  
  prop_data_group <- prop_data[["patient_group"]] %>%
    unique() %>%
    set_names()
  
  map(prop_data_group, function(g) {
    
    test_data <- prop_data %>%
      mutate(test_group = ifelse(patient_group == g,
                                 "target", "rest"))
    
    wilcox.test(corr_p_value ~ test_group, 
                data = test_data,
                alternative = "less") %>%
      broom::tidy()
  }) %>% enframe("patient_group") %>%
    unnest()
  
}

spark_single_slide_ora <- read_csv("./results/spark_ora/spark_ora.csv") %>%
  left_join(sample_dict) %>%
  dplyr::select(gset, corr_p_value, patient_group, sample_id, patient_id) %>%
  mutate(corr_p_value = -log10(corr_p_value)) %>%
  mutate(corr_p_value = ifelse(is.infinite(corr_p_value), 
                               max(corr_p_value),corr_p_value))

spark_ora_wilcox <- spark_single_slide_ora %>%
  dplyr::select(-c("sample_id", "patient_id")) %>%
  group_by(gset) %>% 
  nest() %>%
  mutate(wilcox_res = map(data, run_wilcox_all)) %>%
  dplyr::select(gset, wilcox_res) %>%
  unnest() %>%
  ungroup() %>%
  group_by(patient_group) %>%
  mutate(corr_pval = p.adjust(p.value))
  
filtered_gsets <- spark_ora_wilcox %>%
  dplyr::filter(p.value  <= 0.05) %>%
  arrange(patient_group)

muscle_sets <- filtered_gsets %>%
  dplyr::filter(map_lgl(gset, grepl, pattern = "muscle", ignore.case = T) |
                  map_lgl(gset, grepl, pattern = "contraction", ignore.case = T)  )

ecm_sets <- filtered_gsets %>%
  dplyr::filter(map_lgl(gset, grepl, pattern = "matrisome", ignore.case = T) |
                  map_lgl(gset, grepl, pattern = "ecm", ignore.case = T) | 
                  map_lgl(gset, grepl, pattern = "collagen", ignore.case = T) )

immune_sets <- filtered_gsets %>%
  dplyr::filter(map_lgl(gset, grepl, pattern = "immune", ignore.case = T) |
                  map_lgl(gset, grepl, pattern = "antigen", ignore.case = T) | 
                  map_lgl(gset, grepl, pattern = "neutrophil", ignore.case = T) )

death_sets <- c("REACTOME_PROGRAMMED_CELL_DEATH")

plot_sets <- c(pull(muscle_sets, gset),
               pull(immune_sets, gset),
               death_sets,
               pull(ecm_sets, gset)) %>%
  unique()


plot_dat <- spark_single_slide_ora %>%
  dplyr::filter(gset %in% plot_sets) %>%
  dplyr::mutate(gset = factor(gset,
                              levels = plot_sets)) %>%
  dplyr::group_by(gset, patient_id) %>%
  summarise(corr_p_value = max(corr_p_value)) %>%
  dplyr::mutate(patient_id = factor(patient_id,
                                    levels = pat_order))

write_csv(plot_dat, file = "./results/sample_comparison/gene_patterns/spark_ora_results.csv")


pdf("./results/sample_comparison/gene_patterns/spark_ora_results.pdf", width = 18, height = 8)

plot_dat_plt <- plot_dat %>%
  mutate(corr_p_value = ifelse(corr_p_value>10, 10,
                               corr_p_value)) %>%
  ggplot(aes(x = patient_id, 
             y = gset, 
             fill = corr_p_value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow") +
  ylab("") +
  xlab("")

plot(plot_dat_plt)

dev.off()

# Summarize by zone

naming_groups <- tibble(patient_group = c("group_1",
                                          "group_2",
                                          "group_3"),
                        major_labl = c("CTRL","IZ","FZ"))

group_order <- c("CTRL","RZ","BZ","IZ","FZ")

# Patient description by zones:
zone_description <- sample_dict %>%
  left_join(naming_groups) %>%
  dplyr::mutate(major_labl = ifelse(grepl("RZ", patient_id), 
                                    "RZ",
                                    major_labl)) %>%
  dplyr::mutate(major_labl = ifelse(grepl("BZ", patient_id), 
                                    "BZ",
                                    major_labl)) %>%
  dplyr::mutate(major_labl = ifelse(sample_id=="AKK003_157775", 
                                    "IZ",
                                    major_labl)) %>%
  dplyr::mutate(major_labl = ifelse(sample_id=="AKK001_157785", 
                                    "FZ",
                                    major_labl)) %>%
  dplyr::select(sample_id, patient_id, patient_group, major_labl)


write_csv(zone_description, file = "./markers/zone_description_visium.csv")



summ_dat <- spark_single_slide_ora %>%
  dplyr::select(-c("patient_group", "patient_id")) %>%
  left_join(zone_description, by = c("sample_id")) %>%
  dplyr::filter(gset %in% plot_sets) %>%
  dplyr::mutate(gset = factor(gset,
                              levels = plot_sets)) %>%
  
  group_by(gset, major_labl) %>%
  summarise(median_p = median(corr_p_value)) %>%
  mutate(major_labl = factor(major_labl,
                             levels = group_order))


summ_plt <- summ_dat %>%
  mutate(median_p = ifelse(median_p>15, 15,
                           median_p)) %>%
  ggplot(aes(x = major_labl, 
             y = gset, 
             fill = median_p)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_fill_gradient2(low = "black", mid = "black", high = "yellow") +
  ylab("") +
  xlab("") +
  coord_equal()

write_csv(summ_dat, "./results/sample_comparison/gene_patterns/summ_spark_ora_results.csv")

pdf("./results/sample_comparison/gene_patterns/summ_spark_ora_results.pdf", width = 15, height = 8)

plot(summ_plt)

dev.off()






















