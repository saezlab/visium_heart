# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generation of immune cell response scores
#' 
#' 1. Calculates scores, adds them as a new assay and plots them
#' 
#' Follows structure of slide_scell_map.R
#' 
#' 

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(cowplot)
library(clustree)
library(xlsx)
library(genesorteR)

source("./visium_exploratory/slide_processing.R")
gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")


# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

# Re-clustering + Wilcoxon
for(slide in slides_ids){
  
  print(slide)
  
  # File definition
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  out_immune = sprintf("./results/single_slide/%s/%s_immune_scores.pdf",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  DefaultAssay(visium_slide) = "SCT"
  
  visium_slide[['Immune']] = CreateAssayObject(counts = get_ImmuneScores(seurat_obj = visium_slide,
                                                                      ImmuneSets,ImmuneMarkers))

  DefaultAssay(visium_slide) = "Immune"
  
  Immune_plot_A = SpatialFeaturePlot(object = visium_slide,
                                   features = rownames(visium_slide)[1:10],
                                   stroke = 0)
  
  Immune_plot_B = SpatialFeaturePlot(object = visium_slide,
                                     features = rownames(visium_slide)[11:21],
                                     stroke = 0)
  
  pdf(file = out_immune, width = 15, height = 14,onefile = T)
  plot(Immune_plot_A)
  plot(Immune_plot_B)
  dev.off()

}

