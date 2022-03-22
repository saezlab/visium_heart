# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates correlations of MISTy pipelines

library(tidyverse)
library(Seurat)
library(viridis)
source("./analysis/utils/misty_pipeline.R")
source("./analysis/utils/spatial_plots_utils.R")

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
r2_filter <- 10

get_directions <- function(misty_out_folder, out_alias, exclude = "pattern") {
  
  states_cors_all <- collect_correlations(misty_out_folder, r2_filter, sample_dict = sample_dict)
  
  states_cors_unfilt <- states_cors_all %>%
    dplyr::filter(grepl(exclude,feature_b)) %>%
    dplyr::mutate(correlation = ifelse(target == feature_b, NA, correlation))
  
  states_cors <- states_cors_all %>%
    dplyr::filter(!grepl(exclude,feature_b))
  
  plts_g <- states_cors %>%
    group_by(target, feature_b, patient_group) %>%
    summarize(mean_cor = mean(correlation)) %>%
    ggplot(., aes(x = target, y = feature_b, fill = mean_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    xlab("") +
    ylab("")
  
  plts <- states_cors %>%
    group_by(target, feature_b) %>%
    summarize(mean_cor = mean(correlation)) %>%
    ggplot(., aes(x = target, y = feature_b, fill = mean_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    xlab("") +
    ylab("")
  
  plts_unfilt_g <- states_cors_unfilt %>%
    group_by(target, feature_b, patient_group) %>%
    summarize(mean_cor = mean(correlation)) %>%
    ggplot(., aes(x = target, y = feature_b, fill = mean_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    xlab("") +
    ylab("") 
  
  plts_unfilt <- states_cors_unfilt %>%
    group_by(target, feature_b) %>%
    summarize(mean_cor = mean(correlation)) %>%
    ggplot(., aes(x = target, y = feature_b, fill = mean_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    xlab("") +
    ylab("") 
  
  pdf(paste0("./results/spatial_cors/", out_alias, "_cors.pdf"), height = 4, width = 6)
  
  plot(plts + coord_equal())
  plot(plts_g + facet_wrap(.~patient_group) + coord_equal())
  plot(plts_unfilt + coord_equal())
  plot(plts_unfilt_g + facet_wrap(.~patient_group) + coord_equal())
  
  write_csv(states_cors_all, paste0("./results/spatial_cors/", out_alias, "_cors.csv"))
  
  dev.off()
  
  states_comp <-states_cors %>%
    group_by(target, feature_b) %>%
    nest() %>%
    mutate(pw_res =  map(data, function(dat) {
      
      ggpubr::compare_means(correlation ~ patient_group,  
                            data = dat,
                            method = "wilcox.test", 
                            alternative = "two.sided")
      
    })) %>%
    select(pw_res) %>%
    unnest()
  
  write_csv(states_comp,  paste0("./results/spatial_cors/", out_alias, "_comps.csv"))
  
  states_comp <- plot_box_pw_cor(cors = states_cors, cor_comps = states_comp, area = F)
  
  states_comp <- states_comp %>%
    unnest(pw_data) %>%
    dplyr::filter(p.adj < 0.1) %>%
    select(target,feature_b, bxplt) %>%
    unique() %>%
    mutate(key = paste0(target,"_",feature_b)) %>%
    ungroup() %>%
    select(key, bxplt)
  
  pdf(paste0("./results/spatial_cors/", out_alias, "_comps.pdf"), height = 5, width = 5)
  
  walk2(states_comp$key, states_comp$bxplt, function(k, plt) {
    
    plot(plt + ggtitle(k))
    
    
  })
  
  dev.off()
  
}


#CM
misty_out_folder <- "./results/state_structure/Stressed_CM"
get_directions(misty_out_folder = misty_out_folder,out_alias = "cm_ct",exclude = "CM.")

# Endo
misty_out_folder <- "./results/state_structure/Endo_ct"
get_directions(misty_out_folder = misty_out_folder,out_alias = "endo_ct",exclude = "Endo.")

# Fib
misty_out_folder <- "./results/state_structure/Fib_Myeloid"
get_directions(misty_out_folder = misty_out_folder,out_alias = "fib_states",exclude = "Fib.")

# Myeloid
misty_out_folder <- "./results/state_structure/Myeloid_ct"
get_directions(misty_out_folder = misty_out_folder,out_alias = "myeloid_ct",exclude = "Myeloid.")
