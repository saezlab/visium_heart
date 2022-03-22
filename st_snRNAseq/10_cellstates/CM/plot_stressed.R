# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates figures of examples of main figures

library(tidyverse)
library(Seurat)
library(viridis)
source("./analysis/utils/misty_pipeline.R")
source("./analysis/utils/spatial_plots_utils.R")

# First do the plots of the slides in Figure

# Main

slide_file <- "AKK002_157781.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

  

cell_states <- c("CM-damaged-CM", "CM-healthy-CM", "CM-intermediate-CM")

walk(cell_states, function(ct) {
  
  print(ct)
  
  DefaultAssay(slide) <- "cell_states_pos"
  
  state_plot <- SpatialFeaturePlot(slide, 
                                   features = ct,
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 1.5) + 
    scale_fill_viridis(option = "D")
  
  pdf(paste0("./results/stressed_CMs/visium_examples/",slide_id,"_",ct,".pdf"), height = 4, width = 4)
  
  plot(state_plot)
  
  dev.off()
  
  
})

# Supplement

slide_file <- "AKK003_157777.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

cell_states <- c("CM-damaged-CM", "CM-healthy-CM", "CM-intermediate-CM")

walk(cell_states, function(ct) {
  
  print(ct)
  
  DefaultAssay(slide) <- "cell_states_pos"
  
  state_plot <- SpatialFeaturePlot(slide, 
                                   features = ct,
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 1.5) + 
    scale_fill_viridis(option = "D")
  
  pdf(paste0("./results/stressed_CMs/visium_examples/",slide_id,"_",ct,".pdf"), height = 4, width = 4)
  
  plot(state_plot)
  
  dev.off()
  
  
})

# Now specific interactions

plot_df <- readRDS("./results/stressed_CMs/MISTy_summary.rds")


# Get comparisons between patient groups that are meaningful
plot_df %>%
  dplyr::select(pw_patgroup) %>%
  unnest() %>%
  dplyr::filter(p.adj < 0.15)

# Now generate examples

# vSMCs

samples <- plot_df %>%
  select(data) %>%
  unnest() %>%
  dplyr::filter(Predictor == "vSMCs") %>%
  arrange(-Importance) %>%
  dplyr::slice(1:5) %>%
  pull(sample) %>%
  unique()

walk(samples, function(s) {
  
  print(s)
  
  # vSMCs
  slide_id <- s
  slide_file <- paste0(s, ".rds")
  state_origin <- "cell_states"
  print(slide_file)
  slide_files_folder <- "./processed_visium/objects/"
  
  # Read spatial transcriptomics data and transform states to be useful for modeling
  slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
    positive_states(., assay = state_origin) %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.1)
  
  DefaultAssay(slide) <- "cell_states_pos"
  
  state_plot <- SpatialFeaturePlot(slide, 
                                   features = "CM-damaged-CM",
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 1.5) +
    scale_fill_viridis(option = "D")
  
  DefaultAssay(slide) <- "c2l"
  
  vsmc_plot <- SpatialFeaturePlot(slide, 
                                  features = "vSMCs", 
                                  max.cutoff = "q99", 
                                  stroke = 0,
                                  pt.size.factor = 1.5, 
                                  crop = T) +
    scale_fill_viridis(option = "A")
  
  
  pdf(paste0("./results/stressed_CMs/vsmc_slides/", slide_id,
      "_stressed_allslide.pdf"), height = 4, width = 4)
  
  plot(state_plot)
  
  dev.off()
  
  pdf(paste0("./results/stressed_CMs/vsmc_slides/", slide_id,
      "_vsmcs_allslide.pdf"), height = 4, width = 4)
  
  plot(vsmc_plot)
  
  dev.off()
  
})

# Take specific examples

# vSMCs
slide_file <- "AKK003_157777.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

spots_to_keep <- get_quadrant(visium_slide = slide, row_divisions = 4, col_divisions = 2,coord = c(3,1))

slide <- slide[, spots_to_keep]

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 3.5) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

vsmc_plot <- SpatialFeaturePlot(slide, 
                                features = "vSMCs", 
                                max.cutoff = "q99", 
                                stroke = 0,
                                pt.size.factor = 3.5, 
                                crop = T) +
  scale_fill_viridis(option = "A")


pdf("./results/stressed_CMs/visium_examples/AKK003_157777_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/AKK003_157777_vSMCs.pdf", height = 4, width = 4)

plot(vsmc_plot)

dev.off()

# Fib
samples <- plot_df %>%
  select(data) %>%
  unnest() %>%
  dplyr::filter(Predictor == "Fib") %>%
  arrange(-Importance) %>%
  dplyr::slice(1:5) %>%
  pull(sample) %>%
  unique()

walk(samples, function(s) {
  
  print(s)
  
  # vSMCs
  slide_id <- s
  slide_file <- paste0(s, ".rds")
  state_origin <- "cell_states"
  print(slide_file)
  slide_files_folder <- "./processed_visium/objects/"
  
  # Read spatial transcriptomics data and transform states to be useful for modeling
  slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
    positive_states(., assay = state_origin) %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.1)
  
  DefaultAssay(slide) <- "cell_states_pos"
  
  state_plot <- SpatialFeaturePlot(slide, 
                                   features = "CM-damaged-CM",
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 1.5) +
    scale_fill_viridis(option = "D")
  
  DefaultAssay(slide) <- "c2l"
  
  vsmc_plot <- SpatialFeaturePlot(slide, 
                                  features = "Fib", 
                                  max.cutoff = "q99", 
                                  stroke = 0,
                                  pt.size.factor = 1.5, 
                                  crop = T) +
    scale_fill_viridis(option = "A")
  
  
  pdf(paste0("./results/stressed_CMs/fib_slides/", slide_id,
             "_stressed_allslide.pdf"), height = 4, width = 4)
  
  plot(state_plot)
  
  dev.off()
  
  pdf(paste0("./results/stressed_CMs/fib_slides/", slide_id,
             "_fib_allslide.pdf"), height = 4, width = 4)
  
  plot(vsmc_plot)
  
  dev.off()
  
})

# Generate plot
slide_file <- "Visium_3_CK281.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

spots_to_keep <- get_quadrant(visium_slide = slide, row_divisions = 1, col_divisions = 2,coord = c(1,1))

slide <- slide[, spots_to_keep]

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 2) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

fib_plot <- SpatialFeaturePlot(slide, 
                                features = "Fib", 
                                max.cutoff = "q99", 
                                stroke = 0,
                                pt.size.factor = 2, 
                                crop = T) +
  scale_fill_viridis(option = "A")

pdf("./results/stressed_CMs/visium_examples/Visium_3_CK281_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/Visium_3_CK281_Fib.pdf", height = 4, width = 4)

plot(fib_plot)

dev.off()

# Myeloid
plot_df %>%
  select(data) %>%
  unnest() %>%
  dplyr::filter(Predictor == "Myeloid") %>%
  arrange(-Importance)

slide_file <- "Visium_12_CK290.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

#roi <- get_ct_spots(slide = slide,ct = "CM",filter_prop = 0.1)
#slide <- slide[,roi]

#spots_to_keep <- get_quadrant(visium_slide = slide, row_divisions = 4, col_divisions = 2,coord = c(3,1))
#slide <- slide[, spots_to_keep]

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 2) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

myeloid_plot <- SpatialFeaturePlot(slide, 
                                   features = "Myeloid", 
                                   max.cutoff = "q99", 
                                   stroke = 0,
                                   pt.size.factor = 2, 
                                   crop = T) +
  scale_fill_viridis(option = "A")

pdf("./results/stressed_CMs/visium_examples/Visium_12_CK290_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/Visium_12_CK290_Myeloid.pdf", height = 4, width = 4)

plot(myeloid_plot)

dev.off()

# Remote zones

slide_file <- "Visium_8_CK286.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 2) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

myeloid_plot <- SpatialFeaturePlot(slide, 
                                   features = "Myeloid", 
                                   max.cutoff = "q99", 
                                   stroke = 0,
                                   pt.size.factor = 2, 
                                   crop = T) +
  scale_fill_viridis(option = "A")

pdf("./results/stressed_CMs/visium_examples/Visium_8_CK286_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/Visium_8_CK286_Myeloid.pdf", height = 4, width = 4)

plot(myeloid_plot)

dev.off()

# Next onw

slide_file <- "Visium_5_CK283.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 1.5) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

myeloid_plot <- SpatialFeaturePlot(slide, 
                                   features = "Myeloid", 
                                   max.cutoff = "q99", 
                                   stroke = 0,
                                   pt.size.factor = 1.5, 
                                   crop = T) +
  scale_fill_viridis(option = "A")

pdf("./results/stressed_CMs/visium_examples/Visium_5_CK283_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/Visium_5_CK283_Myeloid.pdf", height = 4, width = 4)

plot(myeloid_plot)

dev.off()



# Adipo

plot_df %>%
  select(data) %>%
  unnest() %>%
  dplyr::filter(Predictor == "Adipo") %>%
  arrange(-Importance)

slide_file <- "AKK002_157781.rds"
state_origin <- "cell_states"
print(slide_file)
slide_files_folder <- "./processed_visium/objects/"

# Read spatial transcriptomics data and transform states to be useful for modeling
slide_id <- gsub("[.]rds", "", slide_file)
slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)

DefaultAssay(slide) <- "cell_states_pos"

state_plot <- SpatialFeaturePlot(slide, 
                                 features = "CM-damaged-CM",
                                 max.cutoff = "q99", 
                                 min.cutoff = "q1",
                                 stroke = 0,
                                 pt.size.factor = 1.5) +
  scale_fill_viridis(option = "D")

DefaultAssay(slide) <- "c2l"

adipo_plot <- SpatialFeaturePlot(slide, 
                                   features = "Adipo", 
                                   max.cutoff = "q99", 
                                   stroke = 0,
                                   pt.size.factor = 1.5, 
                                   crop = T) +
  scale_fill_viridis(option = "A")

pdf("./results/stressed_CMs/visium_examples/AKK002_157781_stressed.pdf", height = 4, width = 4)

plot(state_plot)

dev.off()

pdf("./results/stressed_CMs/visium_examples/AKK002_157781_Adipo.pdf", height = 4, width = 4)

plot(adipo_plot)

dev.off()

