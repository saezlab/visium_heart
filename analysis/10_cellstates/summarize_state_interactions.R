# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We find a consensus structural signature of cell-states
#' 

library(tidyverse)
library(cowplot)
library(mistyR)

misty_out_folder <- "./results/state_structure/fibroblast"

summarize_all_misty = function(misty_out_folder) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  
  misty_res <- collect_results(misty_outs)
  
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

summarize_all_misty("./results/state_structure/fibroblast")
summarize_all_misty("./results/state_structure/cardiomyocyte")
summarize_all_misty("./results/state_structure/endothelial")













