# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot "important" niches to space

library(tidyverse)
library(Seurat)
library(cowplot)

visium_folder = "./processed_visium/objects/"

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample_name = visium_samples)

# Function that plots interesting clusters
plot_all_niches = function(visium_file, sample_name, cluster_list) {
  
  print(sample_name)
  
  visium_slide <- readRDS(visium_file)
  
  cluster_list <- unique(visium_slide$opt_clust_integrated)
  
  Idents(visium_slide) <- "opt_clust_integrated"

  niche_plt <-  SpatialDimPlot(visium_slide,
                 group.by = "opt_clust_integrated",
                 label = TRUE, 
                 label.size = 0,
                 stroke = 0, 
                 label.box = F) + ggtitle(sample_name)
  
  
  all_niches <- map(set_names(cluster_list), function(niche) {
    
    if(niche %in% Idents(visium_slide)) {
      SpatialDimPlot(visium_slide, 
                     cells.highlight = WhichCells(visium_slide, 
                                                  idents = niche),
                     label = F) + ggtitle(niche)
    } else {
      NULL
    }
  })
  
  all_niches <- plot_grid(plotlist = all_niches, 
                          ncol = 4, align = "hv",
                          labels = sample_name)
  
  
  pdf(file = paste0("./results/niche_mapping/mixed/niche_plts_", sample_name,".pdf"), width = 8, height = 9)
  
  plot(niche_plt)
  
  dev.off()
  
  pdf(file = paste0("./results/niche_mapping/sep/niche_plts_", sample_name,".pdf"), width = 18, height = 11)
  
  plot(all_niches)
  
  dev.off()
  
}

all_plts <- pwalk(visium_df, plot_all_niches, cluster_list = niche_filter)













































