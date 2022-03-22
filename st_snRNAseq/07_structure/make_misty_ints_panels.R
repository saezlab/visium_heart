# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# In this script I will generate the examples of cell types and progeny pathways in space

library(Seurat)
library(tidyverse)
library(viridis)
source("./analysis/utils/spatial_plots_utils.R")


# CT-CT

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_20_CK298_mistyassays.rds")
DefaultAssay(visium_slide) <- "c2l"
target_plt <- SpatialFeaturePlot(visium_slide, features = "Endo", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "c2l_juxta"
predictor_plt <- SpatialFeaturePlot(visium_slide, features = "PC", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/tissue_structure/figure2_panels/Visium_20_endotarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_20_pcpred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_9_CK287_mistyassays.rds")
DefaultAssay(visium_slide) <- "c2l"
target_plt <- SpatialFeaturePlot(visium_slide, features = "Endo", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "c2l_juxta"
predictor_plt <- SpatialFeaturePlot(visium_slide, features = "PC", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/tissue_structure/figure2_panels/Visium_9_endotarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_9_pcpred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

# vSMCs Endo

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_3_CK281_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

target_plt <- SpatialFeaturePlot(visium_slide, features = "Endo", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 2.5) +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "c2l_para"
predictor_plt <- SpatialFeaturePlot(visium_slide, features = "vSMCs", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 2.5) +
  scale_fill_viridis(option = "A")

pdf("./results/tissue_structure/figure2_panels/Visium_3_endotarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_3_vsmcspred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

# Path-Path

#plasma -1

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_7_CK285_mistyassays.rds")

DefaultAssay(visium_slide) <- "progeny"

target_plt <- SpatialFeaturePlot(visium_slide, features = "PI3K", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "C")

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "p53", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "C")

pdf("./results/tissue_structure/figure2_panels/Visium_7_pi3ktarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_7_p53pred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

# Added visium 10

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_10_CK288_mistyassays.rds")

DefaultAssay(visium_slide) <- "progeny"

target_plt <- SpatialFeaturePlot(visium_slide, features = "PI3K", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "C")

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "p53", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "C")

pdf("./results/tissue_structure/figure2_panels/Visium_10_pi3ktarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_10_p53pred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()


# CT-Path

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_19_CK297_mistyassays.rds")

DefaultAssay(visium_slide) <- "c2l"

target_plt <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "progeny"

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "Hypoxia", stroke = 1, max.cutoff = "q99", min.cutoff = "q1")  +
  scale_fill_viridis(option = "C")

pdf("./results/tissue_structure/figure2_panels/Visium_19_CMtarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_19_hypoxiapred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

# Visium 15

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_15_CK293_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

target_plt <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 2) +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "progeny"

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "WNT", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 2)  +
  scale_fill_viridis(option = "C")

pdf("./results/tissue_structure/figure2_panels/Visium_15_CMtarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_15_wntpred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

# Visium 20

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_20_CK298_mistyassays.rds")

DefaultAssay(visium_slide) <- "c2l"

target_plt <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "progeny"

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "WNT", stroke = 1, max.cutoff = "q99", min.cutoff = "q1")  +
  scale_fill_viridis(option = "C")

pdf("./results/tissue_structure/figure2_panels/Visium_20_CMtarget.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_20_wntpred.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()


# Visium 13

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_13_CK291_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 2, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

target_plt <- SpatialFeaturePlot(visium_slide, features = "Lymphoid", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3) +
  scale_fill_viridis(option = "A")

target_plt_2 <- SpatialFeaturePlot(visium_slide, features = "Myeloid", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3) +
  scale_fill_viridis(option = "A")

DefaultAssay(visium_slide) <- "progeny"

predictor_plt <- SpatialFeaturePlot(visium_slide, features = "JAK-STAT", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3)  +
  scale_fill_viridis(option = "C")

predictor_plt_2 <- SpatialFeaturePlot(visium_slide, features = "NFkB", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3)  +
  scale_fill_viridis(option = "C")


pdf("./results/tissue_structure/figure2_panels/Visium_13_Lymphoidpred.pdf", height = 5, width = 4)

plot(target_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_13_myeloidpred.pdf", height = 5, width = 4)

plot(target_plt_2)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_13_jakstattarget.pdf", height = 5, width = 4)

plot(predictor_plt)

dev.off()

pdf("./results/tissue_structure/figure2_panels/Visium_13_nfkbtarget.pdf", height = 5, width = 4)

plot(predictor_plt_2)

dev.off()



