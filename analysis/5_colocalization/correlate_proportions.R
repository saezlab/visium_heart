# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate cell2location results from all slides
#' 
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)


# Import path pointers

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()


# First for each slide we will create metavariables that flag the location of a cell-type in a spot

c2l_assay <- "c2l_props"

cell_props <- map(visium_df$visium_file, function(visium_file) {
  
  print(visium_file)
  
  slide <- readRDS(visium_file)
  cell_props <- GetAssayData(slide, assay = c2l_assay) %>% t()
  
  })

names(cell_props) <- paste0("sample", visium_df$sample)

cell_props <- purrr::reduce(cell_props, rbind)

cor_mat <- cor(cell_props,method = "spearman")

pdf("./results/tissue_structure/colocalization/c2l_props_correlation.pdf", height = 7, width = 8)
plot(Heatmap(cor_mat))
dev.off()


