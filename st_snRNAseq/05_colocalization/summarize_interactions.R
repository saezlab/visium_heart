# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates cell maps that are generalizable

library(tidyverse)
library(cowplot)
library(mistyR)
library(Seurat)
source("./analysis/utils/misty_utilities.R")

misty_out_folder <- "./results/tissue_structure/misty/cell_map/"
misty_outs <- list.files(misty_out_folder, full.names = F)
misty_outs <- set_names(misty_outs, gsub("cm_", "", misty_outs) %>%
                          gsub("_c2l", "", .))

# This is the annotation of samples
annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

sample_dict <- read_csv("./markers/visium_patient_anns_revisions.csv") %>%
  left_join(annotation_names) %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched")))

# Create a summary of MISTy results to navigate

misty_res <- collect_results(paste0(misty_out_folder, misty_outs))

# Collect all importances

sample_importances <- misty_res$importances %>% 
  mutate(sample = strsplit(sample, 'cm_') %>%
           map_chr(., ~.x[[2]]) %>%
           gsub("_c2l", "", .)) %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))

# Collect model performance

# R2 distribution of our markers

R2_data <- misty_res$improvements %>%
  dplyr::filter(measure == "multi.R2") %>%
  dplyr::mutate(sample = gsub("_c2l", "", sample) %>%
                  strsplit(.,split = "cm_") %>%
                  map_chr(., ~ last(.x))) %>%
  dplyr::left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  rename("R2" = value)

cell_order <- R2_data %>% 
  group_by(target) %>%
  summarize(med_value = median(R2)) %>%
  arrange(-med_value) %>%
  pull(target)

cells_R2_tile <- ggplot(R2_data, aes(x = factor(target,
                                                levels = cell_order), 
                                     y = sample, fill = R2)) +
  geom_tile() +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(low = "black", high = "yellow")

cells_R2_box <- ggplot(R2_data, aes(x = factor(target,
                                               levels = cell_order), y = R2)) +
  geom_boxplot() +
  geom_point(aes(color = patient_group)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size =12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Explained variance") +
  xlab("")

cells_R2_box_by_group <- ggplot(R2_data, aes(x = patient_group, y = R2)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(.~target)

write_csv(R2_data, "./results/tissue_structure/colocalization/c2l_R2_results.csv")

pdf("./results/tissue_structure/colocalization/c2l_R2_results.pdf", height = 4, width = 6)

plot(cells_R2_box)

dev.off()

pdf("./results/tissue_structure/colocalization/c2l_R2_results_xtra.pdf", height = 7, width = 7)

plot(cells_R2_tile)
plot(cells_R2_box_by_group)

dev.off()


# Define best performers
R2_value = 10

best_performers <- R2_data %>%
  dplyr::select(target, sample, R2) %>%
  dplyr::filter(R2 >= 10) %>%
  dplyr::mutate(best_performer = T)

# Filter importances based on best performance

sample_importances_filt <-  sample_importances %>%
  left_join(best_performers %>% dplyr::select(-R2),
            by = c("Target" = "target", "sample")) %>%
  dplyr::filter(!is.na(best_performer))

write_csv(sample_importances_filt, file = "./results/tissue_structure/colocalization/sample_importances_filt.csv")

# Now, we need to identify the interactions per view and cell that are the most repeated (median)

summarized_interactions <- sample_importances_filt %>%
  group_by(view, Predictor, Target) %>%
  summarize(median_importance = median(Importance)) %>%
  ungroup() %>%
  group_by(view)

# Calculate the importance median significance with wilcoxon tests

importance_test <- sample_importances_filt %>%
  na.omit() %>%
  dplyr::select(view, Predictor, Target, Importance) %>%
  group_by(view, Predictor, Target) %>%
  nest() %>%
  mutate(wres = map(data, function(dat) {
    wilcox.test(dat[[1]], y = NULL, mu = 0, alternative = "greater") %>%
      broom::tidy()
  })) %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  group_by(view) %>%
  mutate(p_adjust = p.adjust(p.value)) %>%
  mutate("sign_label" = ifelse(p_adjust <= 0.15, "*", ""))

# Add info to summarized interactions

summarized_interactions <- summarized_interactions %>%
  left_join(importance_test, by = c("view", "Predictor", "Target"))

# This is just data manipulation

get_sample_order <- function(view_importance) {
  
  #Predictors in rows
  view_importance_mat <- view_importance %>%
    dplyr::select(Predictor, Target, median_importance) %>%
    pivot_wider(names_from = Target, values_from = median_importance, values_fill = 0) %>%
    column_to_rownames("Predictor") %>%
    as.matrix()
  
  #Predictor order
  predictor_order <- hclust(dist(view_importance_mat))
  predictor_order <- predictor_order$labels[predictor_order$order]
  
  #Target order
  target_order <- hclust(dist(t(view_importance_mat)))
  target_order <- target_order$labels[target_order$order]
  
  return(list(predictor = predictor_order, target = target_order))
}

plot_importance <- function(view_importance, cell_order, color_fill = "#8DA0CB") {
  
  view_plt = view_importance %>%
    mutate(Predictor = factor(Predictor, levels =  cell_order$target),
           Target = factor(Target, levels = cell_order$target)) %>%
    ggplot(aes(x = Target, 
               y = Predictor, 
               fill = median_importance,
               label = sign_label)) +
    geom_tile() +
    geom_text() +
    scale_fill_gradient2(high = color_fill, 
                         midpoint = 0,
                         low = "white",
                         na.value = "grey") +
    ggplot2::coord_equal() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size =12),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  return(view_plt)
}

# Per view generate a heatmap

summarized_interactions <- summarized_interactions %>%
  nest()

cell_order <- get_sample_order((summarized_interactions %>%
  dplyr::filter(view == "intra") %>%
  pull(data))[[1]])

intra_importances <- plot_importance((summarized_interactions %>%
                   dplyr::filter(view == "intra") %>%
                   pull(data))[[1]], 
                cell_order = cell_order,
                color_fill = "darkgreen")

juxta_importances <- plot_importance((summarized_interactions %>%
                   dplyr::filter(view == "juxta_5") %>%
                   pull(data))[[1]], 
                cell_order = cell_order,
                color_fill = "orange")

para_importances <- plot_importance((summarized_interactions %>%
                   dplyr::filter(view == "para_15") %>%
                   pull(data))[[1]], 
                cell_order = cell_order)

pdf("./results/tissue_structure/colocalization/misty_importances_ct_intra.pdf", height = 5, width = 5)

plot(intra_importances)

dev.off()

pdf("./results/tissue_structure/colocalization/misty_importances_ct_juxta.pdf", height = 5, width = 5)

plot(juxta_importances)

dev.off()

pdf("./results/tissue_structure/colocalization/misty_importances_ct_para.pdf", height = 5, width = 5)

plot(para_importances)

dev.off()

summarized_interactions %>% unnest() %>% write_csv(., file = "./results/tissue_structure/colocalization/misty_importances_ct.csv")

# What if we use the mean?

pdf("./results/tissue_structure/colocalization/misty_modelperf_ct.pdf")

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

dev.off()

# Supplements -------------------------------------

# Let's separate importances by different conditions and views

summarized_interactions_group <- sample_importances_filt %>%
  group_by(view, Predictor, Target, patient_group) %>%
  summarize(median_importance = median(Importance)) %>%
  ungroup() %>%
  group_by(view)

write_csv(summarized_interactions_group, file = "./results/tissue_structure/colocalization/misty_importances_ct_bygroups.csv")

summarized_interactions_group_plt <- summarized_interactions_group %>%
  ggplot(aes(x = Target, y = Predictor, fill = median_importance)) +
  geom_tile() +
  scale_fill_gradient2(high = "blue", 
                       midpoint = 0,
                       low = "white",
                       na.value = "grey") +
  #ggplot2::coord_equal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(view ~ patient_group, scales = "free") 

pdf("./results/tissue_structure/colocalization/misty_importances_ct_bygroups.pdf")

plot(summarized_interactions_group_plt)

dev.off()
