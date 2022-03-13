# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Evaluates which cell-types are the most important explaining
#' the damaged CM phenotype
library(tidyverse)
library(Seurat)
library(mistyR)
library(ggpubr)
library(viridis)
source("./analysis/utils/misty_utilities.R")
source("./analysis/utils/misty_pipeline.R")
source("./analysis/utils/spatial_plots_utils.R")

# Main ------------------------------------------------------------------------
# Getting sample annotations --------------------------------------------------
sample_dict <- read_csv("./markers/visium_patient_anns_revisions.csv")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Targets
# This is wmean in ROI of interest
run_state_ppline_ct_simple(ROI_ct = "CM",
                           ROI_prop = 0.1,
                           mask_by_prop = F,
                           mask_threshold = 0.1,
                           folder_label = "Stressed_CM",
                           targets = c("CM-damaged-CM"),
                           target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/Stressed_CM/"

misty_out <- mistyR::collect_results(list.dirs(misty_out_folder,full.names = T, recursive = F))

# Which slides are best performers (10 %)

R2_model <- misty_out$improvements %>%
  dplyr::filter(measure == "multi.R2") %>% 
  dplyr::mutate(sample = strsplit(sample, "/") %>%
                  map_chr(., ~ last(.x))) %>%
  dplyr::mutate(sample = gsub("mstate_", "", sample)) %>%
  dplyr::filter(value >= 10) %>%
  arrange(-value)

best_perf <- R2_model %>% pull(sample)

# Check the importances of both views only for best performers

all_importances <- misty_out$importances %>%
  dplyr::mutate(sample = strsplit(sample, "/") %>%
                  map_chr(., ~ last(.x))) %>%
  dplyr::mutate(sample = gsub("mstate_", "", sample)) %>%
  dplyr::filter(sample %in% best_perf,
                view != "intra")

# Get the order based on the median per predictor in all views

ct_order <- all_importances %>%
  group_by(Predictor) %>%
  summarize(med_importance = median(Importance)) %>%
  arrange(-med_importance) %>%
  pull(Predictor)

all_importances_plt <- all_importances %>%
  mutate(Predictor = factor(Predictor, levels = ct_order),
         view = ifelse(view == "intra_pred", "colocalization", "local neighborhood")) %>%
  ggplot(aes(x = Predictor, y = Importance)) +
    geom_boxplot() +
    geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  facet_wrap(view ~ ., ncol = 1)

pdf("./results/stressed_CMs/sup_figures/misty_model_imps.pdf",height = 4, width = 4)

plot(all_importances_plt)

all_importances %>%
  mutate(Predictor = factor(Predictor, levels = ct_order),
         view = ifelse(view == "intra_pred", "colocalization", "local neighborhood")) %>%
  write_csv("./results/stressed_CMs/sup_figures/misty_model_imps.csv")

dev.off()


# Heatmaps

# General
importance_hmap <- all_importances %>%
  mutate(Predictor = factor(Predictor, levels = ct_order),
         view = ifelse(view == "intra_pred", "colocalization", "local neighborhood")) %>%
  group_by(view, Predictor) %>%
  summarize(mean_imp = mean(Importance)) %>%
  ggplot(aes(x = Predictor, y = view, fill = mean_imp)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  coord_equal() +
  scale_fill_gradient2() +
  ylab("") +
  xlab("")

pdf("./results/stressed_CMs/sup_figures/misty_model_imps_hmap.pdf",height = 4, width = 4)

plot(importance_hmap)

dev.off()

# By group
importance_hmap_byg <- all_importances %>%
  left_join(sample_dict, by  = c("sample" = "sample_id")) %>%
  mutate(Predictor = factor(Predictor, levels = ct_order),
         view = ifelse(view == "intra_pred", "colocalization", "local neighborhood")) %>%
  group_by(view, Predictor, patient_group) %>%
  summarize(mean_imp = mean(Importance)) %>%
  ggplot(aes(x = Predictor, y = patient_group, fill = mean_imp)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  scale_fill_gradient2() +
  facet_wrap(. ~ view,ncol = 1) +
  ylab("") +
  xlab("")


pdf("./results/stressed_CMs/sup_figures/misty_model_imps_hmap_byg.pdf",height = 4, width = 4)

all_importances %>%
  left_join(sample_dict, by  = c("sample" = "sample_id")) %>%
  mutate(Predictor = factor(Predictor, levels = ct_order),
         view = ifelse(view == "intra_pred", "colocalization", "local neighborhood")) %>%
  group_by(view, Predictor, patient_group) %>%
  summarize(mean_imp = mean(Importance)) %>%
  write_csv("./results/stressed_CMs/sup_figures/misty_model_imps_hmap_byg.csv")

plot(importance_hmap_byg)

dev.off()

# Now, per view we want to perform pairwise comparisons between groups and between areas

lab_order <- c("CTRL", "RZ")

pw_tests <- all_importances %>%
  left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  group_by(view, Predictor) %>%
  nest() %>%
  dplyr::mutate(pw_patgroup = map(data, function(dat) {
    
    ggpubr::compare_means(Importance ~ patient_group,  
                  data = dat,
                  method = "wilcox.test", 
                  alternative = "two.sided") %>%
      ungroup() 
  })) %>%
  dplyr::mutate(pw_area = map(data, function(dat) {
    
    dat <- dat %>%
      dplyr::filter(major_labl %in% lab_order) %>%
      dplyr::mutate(major_labl = factor(major_labl, levels = lab_order))
      
    
    ggpubr::compare_means(Importance ~ major_labl,  
                          data = dat,
                          method = "t.test", 
                          alternative = "two.sided")
  })) %>%
  dplyr::select(-data)

# Generate panels of patient groups

pw_tests %>% 
  unnest(pw_patgroup) %>%
  dplyr::filter((group1 == "group_2" ) |
                (group2 == "group_2" )) %>%
  arrange(p) %>%
  dplyr::filter(p.adj < 0.15)

pw_tests %>% 
  unnest(pw_area) %>%
  dplyr::arrange(p)


plot_df <- all_importances %>%
  left_join(sample_dict, by = c("sample" = "sample_id")) %>%
  #dplyr::mutate(major_labl = factor(major_labl, levels = lab_order))  %>%
  group_by(view, Predictor) %>%
  nest() %>%
  left_join(pw_tests) %>%
  mutate(view = ifelse(view == "intra_pred", "colocalization", "local_neighborhood"))
  
plt_bplots <- function(view, Predictor, data, pw_patgroup, pw_area) {
  
  #First do boxplots per pat group 
  
  max_val <- max(data$Importance) + 0.05
  
  group_boxplot_df <- ggplot(data, aes(x = patient_group, y = Importance, color = patient_group)) +
        geom_boxplot() +
        geom_point() +
        ggpubr::stat_pvalue_manual(pw_patgroup, label = "p.adj", 
                                   y.position = max_val, 
                                   step.increase = 0.1,
                                   tip.length = 0.01,size = 3) +
        theme_classic() +
        theme(axis.text = element_text(size = 11),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.5, 
                                         size = 11),
              legend.position = "none") +
        ylab("Importance") +
    ggtitle(paste0(view, "-", Predictor))
    
  
  pdf(paste0("./results/stressed_CMs/pw_misty_group/", view, "_", Predictor, ".pdf"), 
      height = 3.5, width = 2)
  
  plot(group_boxplot_df)
  
  dev.off()
  
  #Then do boxplots per area
  
  data <- data %>%
    dplyr::filter(major_labl %in% lab_order) %>%
    dplyr::mutate(major_labl = factor(major_labl, levels = lab_order))
  
  area_boxplot_df <- ggplot(data, aes(x = major_labl, y = Importance)) +
    geom_boxplot() +
    geom_point() +
    ggpubr::stat_pvalue_manual(pw_area, label = "p.adj", 
                               y.position = max_val, 
                               step.increase = 0.1,
                               tip.length = 0.01,size = 3) +
    theme_classic() +
    theme(axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5, 
                                     size = 11),
          legend.position = "none") +
    ylab("Importance") +
    ggtitle(paste0(view, "-", Predictor))
  
  
  pdf(paste0("./results/stressed_CMs/pw_misty_area/", view, "_", Predictor, ".pdf"), 
      height = 3.5, width = 2)
  
  plot(area_boxplot_df)
  
  dev.off()
  
}

pmap(plot_df, plt_bplots)

saveRDS(plot_df, "./results/stressed_CMs/MISTy_summary.rds")

plot_df %>%
  dplyr::filter(Predictor == "vSMCs", view == "colocalization") %>%
  unnest(data) %>%
  pull(patient_group) %>%
  table()
