# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create guide maps to quantify where cell2location is struggling
#' 

library(Seurat)
library(tidyverse)
library(cowplot)

slides_path <- "./visium_results_manuscript/processed_visium_revisions"

param_df <- tibble(slide = list.files(full.names = T, slides_path),
                   sample_name = gsub("[.]rds", "", list.files(full.names = F, slides_path)))


plot_c2lmaps <- function(slide, sample_name) {
  
  visium_slide <- readRDS(slide)
  visium_slide$c2l_sum_score <- colSums(GetAssayData(visium_slide, assay = "c2l")) 
  v_p <- SpatialFeaturePlot(visium_slide, 
                     features = "c2l_sum_score")
  return(v_p)
  
}

param_df <- param_df %>%
  dplyr::mutate(c2l_maps = map2(slide, sample_name, plot_c2lmaps))

panel <- plot_grid(param_df$c2l_maps[[1]], param_df$c2l_maps[[2]], param_df$c2l_maps[[3]],
          param_df$c2l_maps[[4]], param_df$c2l_maps[[5]], param_df$c2l_maps[[6]],
          param_df$c2l_maps[[7]], param_df$c2l_maps[[8]], ncol = 3)


jpeg(file = "./visium_results_manuscript/spot_plots/c2lmaps.jpeg",
     width = 2000, height = 1500)

print(panel)

dev.off()
