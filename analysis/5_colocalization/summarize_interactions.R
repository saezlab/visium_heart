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

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")

# Here we will extract the importances of all views in all samples
sample_importances <- map(misty_outs, function(x) {
  misty_res <- collect_results(paste0(misty_out_folder, x))[["importances.aggregated"]]
})

sample_importances <- sample_importances %>% 
  enframe() %>% 
  unnest() %>%
  #na.omit() %>%
  left_join(sample_dict, by = c("name" = "sample_id"))

# Now, we need to identify the interactions per view and cell that are the most repeated (median)

summarized_interactions <- sample_importances %>%
  group_by(view, Predictor, Target) %>%
  summarize(median_importance = median(Importance)) %>%
  ungroup() %>%
  group_by(view)

# This is just data manipulation

get_sample_order <- function(view_importance) {
  
  #Predictors in rows
  view_importance_mat <- view_importance %>%
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
    ggplot(aes(x = Target, y = Predictor, fill = median_importance)) +
    geom_tile() +
    scale_fill_gradient2(high = color_fill, 
                         midpoint = 0,
                         low = "white",
                         na.value = "grey") +
    ggplot2::coord_equal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
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
# Panel

importance_panel <- cowplot::plot_grid(intra_importances, 
                                       juxta_importances, 
                                       para_importances, 
                                       align = "hv", ncol = 1, 
                                       labels = c("intra", "juxta", "para"))

pdf("./results/tissue_structure/colocalization/misty_importances_ct.pdf", height = 15, width = 5)

plot(importance_panel)

dev.off()

summarized_interactions %>% unnest() %>% write_csv(., file = "./results/tissue_structure/colocalization/misty_importances_ct.csv")


# What if we use the mean?

misty_res <- collect_results(paste0(misty_out_folder, misty_outs))

pdf("./results/tissue_structure/colocalization/misty_modelperf_ct.pdf")

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

dev.off()

# Let's separate importances by different conditions and views

summarized_interactions_group <- sample_importances %>%
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
