# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we show the examples of specific interactions (For paper)

library(tidyverse)
library(ggpubr)
library(Seurat)
library(viridis)
source("./analysis/utils/spatial_plots_utils.R")


annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

sample_importances_filt <- read.csv("./results/sample_comparison/spatial/sample_importances_filt.csv") %>%
  left_join(annotation_names) %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))) 



plot_importance_boxes <- function(sample_importances_filt, view_sel = "intra", cell_pair) {
  
  my_comparisons <- list( c("myogenic-enriched", "ischemic-enriched"), 
                          c("myogenic-enriched", "fibrotic-enriched"), 
                          c("ischemic-enriched", "fibrotic-enriched") )
  
  
  red_data <- sample_importances_filt %>%
    dplyr::select(name, Predictor, 
                  Target, Importance, 
                  patient_group, view) %>%
    dplyr::filter((Predictor == cell_pair[1] &
                     Target == cell_pair[2]) | 
                    (Predictor == cell_pair[2] &
                       Target == cell_pair[1]),
                  view == view_sel) %>%
    group_by(Predictor, Target) %>%
    nest() %>%
    dplyr::mutate(bplot = map(data, function(dat) {
      
      plt <- ggplot(dat,
             aes(x = patient_group, color = patient_group, y = Importance)) +
        geom_boxplot() +
        geom_point() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              legend.position = "none",
              axis.text = element_text(size = 12)) +
        stat_compare_means(comparisons = my_comparisons, size = 3) + # Add pairwise comparisons p-value
        ggtitle(view_sel) 
      

      return(plt)
      
    }))
  
  
  return(red_data)
  
}

# PC to CM colocalization

# PC to  CM

cm_pc <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                               view_sel = "intra",
                               cell_pair = c("CM", "PC")) %>%
  dplyr::filter(Predictor == "PC", Target == "CM")

pdf("./results/sample_comparison/spatial/PC_interactions/PCpred_CMtar_intra_box.pdf", height = 4.5, width = 2)

plot(cm_pc$bplot[[1]])

write_csv(cm_pc$data[[1]], 
          file = "./results/sample_comparison/spatial/PC_interactions/PCpred_CMtar_intra_box.csv")

dev.off()

# Fibrotic example

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_13_CK291_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

cm <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pc <- SpatialFeaturePlot(visium_slide, features = "PC", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/PC_interactions/Visium_13_chronic_CM.pdf", height = 5, width = 4)

plot(cm)

dev.off()

pdf("./results/sample_comparison/spatial/PC_interactions/Visium_13_chronic_PC.pdf", height = 5, width = 4)

plot(pc)

dev.off()


# Myogenic example

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_5_CK283_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,2))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

cm <- SpatialFeaturePlot(visium_slide, features = "CM", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pc <- SpatialFeaturePlot(visium_slide, features = "PC", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/PC_interactions/Visium_5_CM.pdf", height = 5, width = 4)

plot(cm)

dev.off()

pdf("./results/sample_comparison/spatial/PC_interactions/Visium_5_PC.pdf", height = 5, width = 4)

plot(pc)

dev.off()

# Fibrotic example

# Fib to vSMCs

vsmcs_fib <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                                   view_sel = "juxta_5",
                                   cell_pair = c("Fib", "vSMCs")) %>%
  dplyr::filter(Predictor == "vSMCs", Target == "Fib")


pdf("./results/sample_comparison/spatial/Fib_interactions/vSMCspred_Fibtar_juxta_box.pdf", height = 4.5, width = 2)

plot(vsmcs_fib$bplot[[1]])

write_csv(vsmcs_fib$data[[1]], 
          file = "./results/sample_comparison/spatial/Fib_interactions/vSMCspred_Fibtar_juxta_box.csv")

dev.off()

# Myogenic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_3_CK281_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 1.8) +
  scale_fill_viridis(option = "A")

vsmcs <- SpatialFeaturePlot(visium_slide, features = "vSMCs", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 1.8) +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_3_fib.pdf", height = 5, width = 4)

plot(fib)

dev.off()

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_3_vsmcs.pdf", height = 5, width = 4)

plot(vsmcs)

dev.off()

# Myogenic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_8_CK286_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 3,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 2.5) +
  scale_fill_viridis(option = "A")

vsmcs <- SpatialFeaturePlot(visium_slide, features = "vSMCs", stroke = 1, max.cutoff = "q99", min.cutoff = "q1",pt.size.factor = 2.5) +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_8_fib.pdf", height = 5, width = 4)

plot(fib)

dev.off()

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_8_vsmcs.pdf", height = 5, width = 4)

plot(vsmcs)

dev.off()

# Fibrotic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/AKK004_157772_mistyassays.rds")

DefaultAssay(visium_slide) <- "c2l"

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

vsmcs <- SpatialFeaturePlot(visium_slide, features = "vSMCs", stroke = 1, max.cutoff = "q99", min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/Fib_interactions/157772_fib.pdf", height = 5, width = 4)

plot(fib)

dev.off()

pdf("./results/sample_comparison/spatial/Fib_interactions/157772_vsmcs.pdf", height = 5, width = 4)

plot(vsmcs)

dev.off()

# Ischemic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_19_CK297_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 2, col_divisions = 2,coord = c(2,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

fib <- SpatialFeaturePlot(visium_slide, features = "Fib", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3) +
  scale_fill_viridis(option = "A")

vsmcs <- SpatialFeaturePlot(visium_slide, features = "vSMCs", stroke = 1, max.cutoff = "q99", min.cutoff = "q1", pt.size.factor = 3) +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_19_fib.pdf", height = 5, width = 4)

plot(fib)

dev.off()

pdf("./results/sample_comparison/spatial/Fib_interactions/Visium_19_vsmcs.pdf", height = 5, width = 4)

plot(vsmcs)

dev.off()

# Immune example

# Myeloid to Lymphoid

myel_lymph <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                                   view_sel = "juxta_5",
                                   cell_pair = c("Myeloid", "Lymphoid")) %>%
  dplyr::filter(Predictor == "Lymphoid", Target == "Myeloid")


pdf("./results/sample_comparison/spatial/immune_interactions/lymphpred_Myeltar_juxta_box.pdf", height = 4.5, width = 2)

plot(myel_lymph$bplot[[1]])

write_csv(myel_lymph$data[[1]], 
          file = "./results/sample_comparison/spatial/immune_interactions/lymphpred_Myeltar_juxta_box.csv")

dev.off()

# Ischemic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/AKK002_157779_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 2, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

lymph <- SpatialFeaturePlot(visium_slide, features = "Lymphoid", 
                            stroke = 1, max.cutoff = "q99", 
                            min.cutoff = "q1",pt.size.factor = 3.3) +
  scale_fill_viridis(option = "A")

myel <- SpatialFeaturePlot(visium_slide, features = "Myeloid", 
                           stroke = 1, max.cutoff = "q99", 
                           min.cutoff = "q1",pt.size.factor = 3.3) +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/immune_interactions/AKK002_157779_lymphoid.pdf", height = 5, width = 4)

plot(lymph)

dev.off()

pdf("./results/sample_comparison/spatial/immune_interactions/AKK002_157779_myeloid.pdf", height = 5, width = 4)

plot(myel) 

dev.off()

# Ischemic - 2

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_9_CK287_mistyassays.rds")

#spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

#visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

lymph <- SpatialFeaturePlot(visium_slide, features = "Lymphoid", 
                            stroke = 1, max.cutoff = "q99", 
                            min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

myel <- SpatialFeaturePlot(visium_slide, features = "Myeloid", 
                           stroke = 1, max.cutoff = "q99", 
                           min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/immune_interactions/Visium_9_lymphoid.pdf", height = 5, width = 4)

plot(lymph)

dev.off()

pdf("./results/sample_comparison/spatial/immune_interactions/Visium_9_myeloid.pdf", height = 5, width = 4)

plot(myel) 

dev.off()

# Myogenic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/Visium_4_CK282_mistyassays.rds")

DefaultAssay(visium_slide) <- "c2l"

lymph <- SpatialFeaturePlot(visium_slide, features = "Lymphoid", 
                            stroke = 1, max.cutoff = "q99", 
                            min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

myel <- SpatialFeaturePlot(visium_slide, features = "Myeloid", 
                           stroke = 1, max.cutoff = "q99", 
                           min.cutoff = "q1") +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/immune_interactions/Visium_4_lymphoid.pdf", height = 5, width = 4)

plot(lymph)

dev.off()

pdf("./results/sample_comparison/spatial/immune_interactions/Visium_4_myeloid.pdf", height = 5, width = 4)

plot(myel) 

dev.off()

# Myogenic

visium_slide  <- readRDS("./processed_visium/misty_red_objects/AKK006_157771_mistyassays.rds")

spots_to_keep <- get_quadrant(visium_slide = visium_slide, row_divisions = 2, col_divisions = 2,coord = c(1,1))

visium_slide <- visium_slide[, spots_to_keep]

DefaultAssay(visium_slide) <- "c2l"

lymph <- SpatialFeaturePlot(visium_slide, features = "Lymphoid", 
                            stroke = 1, max.cutoff = "q99", 
                            min.cutoff = "q1", pt.size.factor = 3.3) +
  scale_fill_viridis(option = "A")

myel <- SpatialFeaturePlot(visium_slide, features = "Myeloid", 
                           stroke = 1, max.cutoff = "q99", 
                           min.cutoff = "q1", pt.size.factor = 3.3) +
  scale_fill_viridis(option = "A")

pdf("./results/sample_comparison/spatial/immune_interactions/AKK006_157771_lymphoid.pdf", height = 5, width = 4)

plot(lymph)

dev.off()

pdf("./results/sample_comparison/spatial/immune_interactions/AKK006_157771_myeloid.pdf", height = 5, width = 4)

plot(myel) 

dev.off()


