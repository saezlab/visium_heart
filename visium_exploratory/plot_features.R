# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Visualization of genes using Seurat
#' 
#' This script assumes that your directory contains
#' the single slide folder with the slides'
#' rds files in the following structure and is
#' relative to the location where you run this script
#' 
#' "./results/single_slide/%s/%s.rds"
#' 

require(Seurat)
require(ggplot2)

# Annotation of visium and scell ids if needed
sample_dictionary = readRDS("./sample_dictionary.rds")

# Here you define the slide you are interested
# and read it
slide = "157772" #Your selected slide
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)

# Here you define the type of data you want to plot, by default it
# is gene expression

# To check the available assays run:
# Assays(visium_slide)

chosen_assay = "SCT"

DefaultAssay(visium_slide) = chosen_assay

#Here define your gene or feature

feature = "POSTN"

spatial_plot = SpatialFeaturePlot(object = visium_slide,
                        features = feature)

print(spatial_plot)
  












