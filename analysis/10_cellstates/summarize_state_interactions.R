# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We find a consensus structural signature of cell-states
#' 

library(tidyverse)
library(cowplot)
library(mistyR)

patient_annotation <- readRDS("./markers/visium_patient_anns_revisions.rds")
misty_out_folder <- "./results/state_structure/endothelial"

summarize_all_misty = function(misty_out_folder, r2_filter = 0.2) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .))
  
  misty_res <- collect_results(misty_outs)
  
  model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                  grepl("multi.R2", measure)) %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .)) %>%
    left_join(patient_annotation, by = c("sample" = "sample_id"))
  
  mean_multi <- model_performance %>% 
    dplyr::filter(measure == "multi.R2") %>% 
    group_by(sample,patient_group) %>% 
    summarise(meanR2 = mean(value)) %>% 
    arrange(-meanR2)
  
  selected_samples <- dplyr::filter(mean_multi, meanR2 >= r2_filter) %>%
    pull(sample)
  
  selected_files <- paramdf %>%
    dplyr::filter(sample %in% selected_samples)
  
  misty_res <- collect_results(selected_files$misty_out_file)
  
  pdf(file = paste0(misty_out_folder,"/summary_plot_mean.pdf"))
  
  mistyR::plot_improvement_stats(misty_res)
  mistyR::plot_view_contributions(misty_res)
  
  mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res, "intra_cts", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res, "intra_progeny", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res, "para_cts_15", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res, "para_progeny_15", cutoff = 0)
  
  dev.off()
  
}

summarize_all_misty("./results/state_structure/Fib")
summarize_all_misty("./results/state_structure/cardiomyocyte")
summarize_all_misty("./results/state_structure/endothelial")


model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                         grepl("multi.R2", measure)) %>%
  mutate(sample = strsplit(sample, 'mstate_') %>%
           map_chr(., ~.x[[2]]) %>%
           gsub("_c2l","", .)) %>%
  left_join(patient_annotation, by = c("sample" = "sample_id"))

model_performance %>% 
  dplyr::filter(measure == "multi.R2") %>% 
  group_by(sample,patient_group) %>% 
  summarise(meanR2 = mean(value)) %>% 
  arrange(-meanR2)


ggplot(model_performance, aes(x = patient_group, color = patient_group, y = value)) +
  geom_boxplot() +
  facet_wrap(measure~target)
  
  
  ggplot(., aes(x = target, y = sample, fill = value)) +
  geom_tile() +
  facet_wrap(.~measure) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "#ffd89b", high = "#19547b", limits = c(0,1))












