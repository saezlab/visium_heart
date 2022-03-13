# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# In this script I will generate the specific analysis of niche definition

library(Seurat)
library(tidyverse)
library(viridis)
source("./analysis/utils/spatial_plots_utils.R")


niche_cols <- list(niche_1 = "#D51F26",
                 niche_2 = "#272E6A",
                 niche_3 = "#208A42",
                 niche_4 = "#89288F",
                 niche_5 = "#F47D2B",
                 niche_6 = "#FEE500",
                 niche_7 = "#8A9FD1",
                 niche_8 = "#C06CAB",
                 niche_9 = "#D8A767") %>%
  unlist()

visium_slide <- readRDS("./processed_visium/misty_red_objects/Visium_20_CK298_mistyassays.rds")
spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 5, col_divisions = 3,coord = c(5,2))

visium_slide <- visium_slide[,spots_to_keep]

DefaultAssay(visium_slide) <- "c2l_props"

niche_map <- SpatialDimPlot(visium_slide,group.by = "opt_clust_integrated", pt.size.factor = 4) +
  scale_fill_manual(values = niche_cols)

myeloid <- SpatialFeaturePlot(visium_slide, features = "Myeloid", stroke = 1, pt.size.factor = 4) +
  scale_fill_viridis(option = "A") 

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, pt.size.factor = 4) +
  scale_fill_viridis(option = "A")

cm <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, pt.size.factor = 4, max.cutoff = "q99") +
  scale_fill_viridis(option = "A")

histology <- SpatialFeaturePlot(visium_slide, features = "Myeloid",max.cutoff = "q99", stroke = 1, alpha = 0, crop = T, pt.size.factor = 5) +
  scale_fill_viridis(option = "A") 

pdf("./results/niche_mapping/niche_example/niche_map_visium_20.pdf",height = 4, width = 4)

plot(niche_map)

dev.off()

pdf("./results/niche_mapping/niche_example/myeloidprops_visium_20.pdf",height = 4, width = 4)

plot(myeloid)

dev.off()

pdf("./results/niche_mapping/niche_example/fibprops_visium_20.pdf",height = 4, width = 4)

plot(fib)

dev.off()

pdf("./results/niche_mapping/niche_example/cmprops_visium_20.pdf",height = 4, width = 4)

plot(cm)

dev.off()


visium_slide <- readRDS("./processed_visium/misty_red_objects/Visium_18_CK296_mistyassays.rds")
spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 2, col_divisions = 2,coord = c(2,2))

visium_slide <- visium_slide[,spots_to_keep]

DefaultAssay(visium_slide) <- "c2l_props"

niche_map <- SpatialDimPlot(visium_slide,group.by = "opt_clust_integrated", pt.size.factor = 3) +
  scale_fill_manual(values = niche_cols)

myeloid <- SpatialFeaturePlot(visium_slide, features = "Myeloid", stroke = 1, pt.size.factor = 3, max.cutoff = "q99") +
  scale_fill_viridis(option = "A") 

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, pt.size.factor = 3, max.cutoff = "q99") +
  scale_fill_viridis(option = "A")

cm <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, pt.size.factor = 3, max.cutoff = "q99") +
  scale_fill_viridis(option = "A")

histology <- SpatialFeaturePlot(visium_slide, features = "Myeloid",max.cutoff = "q99", stroke = 1, alpha = 0, crop = T, pt.size.factor = 5) +
  scale_fill_viridis(option = "A") 

pdf("./results/niche_mapping/niche_example/niche_map_visium_18.pdf",height = 4, width = 4)

plot(niche_map)

dev.off()

pdf("./results/niche_mapping/niche_example/myeloidprops_visium_18.pdf",height = 4, width = 4)

plot(myeloid)

dev.off()

pdf("./results/niche_mapping/niche_example/fibprops_visium_18.pdf",height = 4, width = 4)

plot(fib)

dev.off()

pdf("./results/niche_mapping/niche_example/cmprops_visium_18.pdf",height = 4, width = 4)

plot(cm)

dev.off()
