# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Test a simplified version of the spatial analysis of interacting cells of interest

library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")
source("./analysis/utils/misty_pipeline.R")

future::plan(future::multisession)

# Get patient annotation
annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

sample_dict <- read_csv("./markers/visium_patient_anns_revisions.csv") %>%
  left_join(annotation_names) %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched")))

# Explained variance filter
r2_filter <- 10
misty_out_folder <- "./results/state_structure/CM_ct_healthy/"
misty_outs <- list.files(misty_out_folder, full.names = T)
misty_outs <- misty_outs[grepl("mstate", misty_outs)]

misty_res <- collect_results(misty_outs)

model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                grepl("multi.R2", measure)) %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]]))

R2_data <- model_performance %>% 
  dplyr::filter(measure == "multi.R2") %>% 
  group_by(target) %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))

# First show the distribution among samples

R2_data_plt <- R2_data %>%
  ggplot(aes(x = patient_group, y = value, color = patient_group)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("Explained variance") +
  xlab("")

R2_data_plt <- R2_data %>%
  ggplot(aes(x = major_labl, y = value, color = major_labl)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("Explained variance") +
  xlab("")

# Make 

best_performers <- R2_data %>% 
  dplyr::filter(value >= r2_filter) %>%
  pull(sample) 

importances_filtered <- misty_res$importances %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]])) %>%
  left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  dplyr::filter(sample %in% best_performers,
                view != "intra") 

intra_order <- importances_filtered %>%
  dplyr::filter(view == "intra_pred") %>%
  group_by(Predictor) %>%
  summarise(median_imp = median(Importance)) %>%
  arrange(-median_imp) %>%
  pull(Predictor)

intra_plt <- importances_filtered %>%
  dplyr::filter(view == "intra_pred") %>%
  ggplot(aes(x = factor(Predictor,
                        levels = intra_order), y = Importance)) +
  geom_boxplot() +
  theme_classic() +
  geom_point() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("")

my_comparisons <- list( c("myogenic-enriched", "ischemic-enriched"), 
                        c("myogenic-enriched", "fibrotic-enriched"), 
                        c("ischemic-enriched", "fibrotic-enriched") )


plts <- map(c("vSMCs","Fib",
              "Adipo", "Myeloid"), function(ct) {
                
                print(ct)
                intra_ct_plt <- importances_filtered %>%
                  dplyr::filter(view == "intra_pred",
                                Predictor %in% ct) %>%
                  ggplot(aes(x = patient_group, color = patient_group, y = Importance)) +
                  geom_boxplot() +
                  theme_classic() +
                  geom_point() +
                  theme(axis.text = element_text(size = 10),
                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                  xlab("") +
                  stat_compare_means(comparisons = my_comparisons, label="p.adj")
              })




cors_filtered <- paramdf %>%
  dplyr::select(sample, cor_res) %>%
  unnest() %>%
  dplyr::rename("correlation" = value) %>%
  dplyr::mutate(feature_a = gsub("-", ".", feature_a),
                feature_b = gsub("-", ".", feature_b)) %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))

















