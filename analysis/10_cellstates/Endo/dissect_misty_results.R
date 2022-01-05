# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Test a simplified version of the spatial analysis of interacting cells of interest

library(tidyverse)
library(Seurat)
library(mistyR)

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
r2_filter <- .10
misty_out_folder <- "./results/state_structure/Endo_ct/"
misty_outs <- list.files(misty_out_folder, full.names = T)
misty_outs <- misty_outs[grepl("mstate", misty_outs)]

misty_res <- collect_results(misty_outs)

model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                grepl("multi.R2", measure) |
                                                                grepl("multi.RMSE", measure))  %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]]))

R2_data <- model_performance %>% 
  dplyr::filter(measure == "multi.R2") %>% 
  group_by(target) %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))

RMSE_data <- model_performance %>% 
  dplyr::filter(measure == "multi.RMSE") %>% 
  group_by(target) %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))


# First show that there are differences in performance of different markers
mrkr_order <- R2_data %>%
  group_by(target) %>%
  summarize(med_imp = median(value)) %>%
  arrange(-med_imp) %>%
  pull(target)

R2_data_plt_mrkrs <- R2_data %>%
  ggplot(aes(x = factor(target,
                        levels = mrkr_order), y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("Explained variance") +
  xlab("") +
  stat_compare_means(label.y = 0.9)

pdf(paste0(misty_out_folder,"mrkr_R2.pdf"), height = 4, width = 3)
plot(R2_data_plt_mrkrs)
dev.off()


RMSE_data_plt_mrkrs <- RMSE_data %>%
  ggplot(aes(x = factor(target,
                        levels = mrkr_order), y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("RMSE") +
  xlab("") +
  stat_compare_means(label.y = 4)

pdf(paste0(misty_out_folder,"mrkr_RMSE.pdf"), height = 4, width = 3)
plot(R2_data_plt_mrkrs)
dev.off()

# Then show the contributions

contributions <- misty_res$contributions %>%
  dplyr::filter(!stringr::str_starts(.data$view, "p\\.") &
                  .data$view != "intercept") %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]])) %>%
  group_by(sample) %>%
  nest() %>%
  dplyr::mutate(data = map(data, function(dat) {
    dat %>%
      dplyr::group_by(.data$target, .data$view) %>%
      dplyr::summarise(mean = mean(.data$value), .groups = "drop_last") %>%
      dplyr::mutate(fraction = abs(.data$mean) / sum(abs(.data$mean))) %>%
      dplyr::ungroup()
  })) %>%
  unnest() %>%
  ungroup()

cont_plt_mrkrs <- ggplot(contributions, aes(x = target, y = fraction, fill = view)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("View contribution") +
  xlab("")

write_csv(contributions, file = paste0(misty_out_folder,"mrkr_contributions.csv"))

pdf(paste0(misty_out_folder,"mrkr_contributions.pdf"), height = 4, width = 5)
plot(cont_plt_mrkrs)
dev.off()

# Then get best in class

best_performers <- R2_data %>% 
  dplyr::filter(value >= r2_filter) %>%
  pull(sample) 

importances_filtered <- misty_res$importances %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]])) %>%
  left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  dplyr::filter(sample %in% best_performers) 

# Plot importances for all the slides, regardless the patient group

summarized_importances <- importances_filtered %>%
  na.omit() %>%
  group_by(view, Target, Predictor) %>%
  summarise(mean_imp = mean(Importance),
            median_imp = median(Importance)) %>%
  ungroup()


# Calculate difference from 0

imp_p <- importances_filtered %>%
  na.omit() %>%
  select(view, Predictor, Target, Importance) %>%
  group_by(view, Predictor, Target) %>%
  nest() %>%
  mutate(wilcox_res = map(data, function(dat) {
    
    wilcox.test(dat$Importance, y = NULL, mu = 0, alternative = "greater") %>%
      broom::tidy()
    
    
  })) %>%
  dplyr::select(-data) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(p_corr = p.adjust(p.value)) %>%
  dplyr::mutate(sign_symbol = ifelse(p_corr <= 0.05, "*", "")) %>%
  dplyr::select(-c("statistic", "method", "alternative"))

summarized_importances <- summarized_importances %>%
  left_join(imp_p, by = c("view", "Predictor", "Target"))

write_csv(summarized_importances, file = paste0(misty_out_folder,"summarized_importances.csv"))

# Plot heatmaps

summarized_importances_plts <- summarized_importances %>%
  group_by(view) %>%
  nest() %>%
  mutate(gplots = map2(view, data, function(v, dat) {
    
    imp_plot <- dat %>%
      ggplot(., aes(x = Target, y = Predictor, fill = median_imp, label = sign_symbol)) +
      geom_tile() +
      geom_text() +
      theme_classic() +
      scale_fill_gradient2(high = "#8DA0CB", 
                           midpoint = 0,
                           low = "white",
                           na.value = "grey") +
      ggplot2::coord_equal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(v)
    
  }))

pdf(file = paste0(misty_out_folder,"/median_importances.pdf"), height = 4, width = 4)

walk(summarized_importances_plts$gplots, plot)

dev.off()

# Now we need to ask about each single marker
# Do we see a difference in R2, contributions?

my_comparisons <- list( c("myogenic-enriched", "ischemic-enriched"), 
                        c("myogenic-enriched", "fibrotic-enriched"), 
                        c("ischemic-enriched", "fibrotic-enriched") )

R2_group_comparisons <- R2_data %>%
  group_by(target) %>%
  nest() %>%
  mutate(gplot = map2(target, data, function(trgt, dat) {
    
    ggplot(dat, aes(x = patient_group, y = value, color = patient_group)) +
      geom_boxplot() +
      geom_point() +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(trgt) +
      stat_compare_means(comparisons = my_comparisons, label="p.adj") +
      ylab("Explained Variance") +
      xlab("")
    
  }))
  

pdf(file = paste0(misty_out_folder,"/R2_group_comparisons.pdf"), height = 5, width = 4)

walk(R2_group_comparisons$gplot, plot)

dev.off()

# Now contributionss
# Space may matter more in one group vs other

contributions_filtered <- contributions %>%
  left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  dplyr::filter(sample %in% best_performers) %>%
  group_by(target, view) %>%
  nest() 

contribution_gcomp_plts <- pmap(contributions_filtered, function(target, view, data) {
    
    ggplot(data, aes(x = patient_group, y = fraction, 
                    color = patient_group)) +
      geom_boxplot() +
      geom_point() +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(target) +
      stat_compare_means(comparisons = my_comparisons, label="p.adj") +
      ylab("View contribution") +
      xlab(view)
    
    
  })

pdf(file = paste0(misty_out_folder,"/contribution_group_comparisons.pdf"), height = 5, width = 4)

walk(contribution_gcomp_plts, plot)

dev.off()

# Finally importances
# Here we are only going to compare importances of capillary cells

target <- "Endo.Capillary.Endo"

view_name <- "para_pred_5"

imp_comp_plts <- importances_filtered %>%
  dplyr::filter(view == view_name,
                Target == target) %>%
  group_by(Predictor) %>%
  nest() %>%
  mutate(gplot = map2(Predictor, data, function(pred, dat) {
    
    ggplot(dat, aes(x = patient_group, y = Importance, color = patient_group)) +
      geom_boxplot() +
      geom_point() +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(pred) +
      stat_compare_means(comparisons = my_comparisons, label="p.adj") +
      ylab("Importance") +
      xlab("")
    
  }))


pdf(file = paste0(misty_out_folder,"/",target, "_", view_name,"imp_comp.pdf"), height = 5, width = 4)

walk(imp_comp_plts$gplot, plot)

dev.off()












