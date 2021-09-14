# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Plot the c2l plots to have a visual representation of the cells

library(Seurat)
library(tidyverse)
library(cowplot)

visium_folder = "./processed_visium/objects/"
out_folder_abundance = "results/deconvolution_models/plots/abundance/"
out_folder_proportion = "results/deconvolution_models/plots/proportion/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

assay_names <- set_names(c("c2l", "c2l_props"))

SpatialPal = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) 

walk2(visium_df$visium_file, visium_df$sample, function(f_path, s) {
  
  print(s)
  
  visium_slide <- readRDS(f_path)
  
  assay_plots <- map(assay_names, function(assay_name){
    
    features <- rownames(GetAssayData(visium_slide,
                                      slot = "data",
                                      assay = assay_name))
    
    DefaultAssay(visium_slide) <- assay_name
    
    feat_plots <- map(features, function(spec_f) {
      
      f_p <- SpatialFeaturePlot(visium_slide, 
                                features = spec_f)
      
      if(grepl(pattern = "props", assay_name)) {
        f_p <- f_p + 
          scale_fill_gradientn(colours = SpatialPal(length(seq(0,1,.1))),
                               limits = c(0,1),
                               breaks = c(0, 0.25, 0.5, 0.75, 1))
      }
      
      return(f_p)
      
    })
    
    niche_plots <- plot_grid(plotlist = feat_plots, ncol = 4)
    
  })
  
  # Plot abundances
  
  abundance_f <- paste0(out_folder_abundance, s,".pdf")
  
  pdf(file = abundance_f, width = 15, height = 15)
  
  plot(assay_plots[[1]])
  
  dev.off()
  
  prop_f <- paste0(out_folder_proportion, s,".pdf")
  
  pdf(file = prop_f, width = 15, height = 15)
  
  plot(assay_plots[[2]])
  
  dev.off()
  
})
