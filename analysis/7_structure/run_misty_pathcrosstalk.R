# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run MISTy for spatial interactions between slides
#' 
library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")

future::plan(future::multisession)

# Pipeline definition:
run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/tissue_structure/misty/pathway_map/pm_") {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay,
                      "intra_ct" = "c2l",
                      "para_ct" = "c2l")
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features,
                        "intra_ct" = NULL,
                        "para_ct" = NULL)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para",
                     "intra_ct" = "intra",
                     "para_ct" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 5,
                      "para" = 15,
                      "intra_ct" = NULL,
                      "para_ct" = 15)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}

# Main ------------------------------------------------------------------------

# Getting parameters

slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Progeny paths
assay_label <- "progeny"

misty_outs <- map(slide_files, function(slide_file){
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  
  mout <- run_colocalization(slide = slide,
                             useful_features = useful_features,
                             out_label = slide_id,
                             assay = assay,
                             misty_out_alias = "./results/tissue_structure/misty/pathway_map/pm_")
  
  misty_res_slide <- collect_results(mout)
  
  plot_folder <- paste0(mout, "/plots")
  
  system(paste0("mkdir ", plot_folder))
  
  pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
  
  mistyR::plot_improvement_stats(misty_res_slide)
  mistyR::plot_view_contributions(misty_res_slide)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "juxta_5", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "juxta_5", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "para_15", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res_slide, "para_15", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "intra_ct", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res_slide, "para_ct_15", cutoff = 0)
  
  dev.off()
  
  return(mout)
  
})

# Now get the summary of this ------------------------------------------------------------------------

misty_outs <- list.dirs("./results/tissue_structure/misty/pathway_map",full.names = T,recursive = F)
misty_res <- collect_results(misty_outs)

pdf("./results/tissue_structure/misty_modelperf_progeny_ct.pdf", height = 6, width = 6)

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)
mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
mistyR::plot_interaction_heatmap(misty_res, "para_15", cutoff = 0)
mistyR::plot_interaction_heatmap(misty_res, "intra_ct", cutoff = 0)
mistyR::plot_interaction_heatmap(misty_res, "para_ct_15", cutoff = 0)

dev.off()

median_importances <- misty_res$importances %>%
  group_by(view, Predictor, Target) %>%
  summarize(median_importance = median(Importance)) %>%
  ungroup() %>%
  group_by(view) %>%
  nest()

median_importances %>% 
  unnest() %>%
  write_csv(., file = "./results/tissue_structure/misty_modelperf_progeny_ct_median.csv")



# Some plotting and manipulation functions

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
  
plot_importance_eqassay <- function(view_importance, cell_order, color_fill = "#8DA0CB") {
  
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

plot_importance_difassay <- function(view_importance, cell_order, color_fill = "#8DA0CB") {
  
  view_plt = view_importance %>%
    mutate(Predictor = factor(Predictor, levels =  cell_order$predictor),
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

# Generate panels

pdf("./results/tissue_structure/misty_modelperf_progeny_ct_median.pdf", height = 6, width = 6)

walk2(median_importances$view, median_importances$data, function(v, dat) {
  
  print(v) 
  
  cell_order <- get_sample_order(dat)
  
  if(grepl("ct", v)) {
    
    plt <- plot_importance_difassay(view_importance = dat,
                             cell_order = cell_order) +
      ggtitle(v)
    
  } else {
    
    plt <- plot_importance_eqassay(view_importance = dat,
                                    cell_order = cell_order) +
      ggtitle(v)
    
  }
  
  plot(plt)
  
} )

dev.off()









