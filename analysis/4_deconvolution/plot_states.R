# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we will visualize the states of each cell-type in spots
#' where the proportion is higher than a human threshold

library(Seurat)
library(tidyverse)
library(cowplot)
source("./analysis/utils/funcomics.R")

# Given a proportion threshold and a cell-type we return a reduced matrix
# same number of spots but filled with NANs where info is not found

mask_state_matrix <- function(cell_props, cell_states, 
                              cell_type, cell_type_prop) {
  
  cell_props[cell_props <= cell_type_prop] = NA
  cell_props[!is.na(cell_props)] = 1
  
  cell_states_ix <- grepl(cell_type, colnames(cell_states))
  cell_states[, cell_states_ix] <- cell_states[, cell_states_ix] * cell_props[, cell_type]
  
  return(cell_states[, cell_states_ix])
  
}

filter_state_assay <- function(cts, 
                               slide, 
                               ref_assay = "c2l_props",
                               state_assay = "cell_states") {
  
  cell_props <- GetAssayData(slide, assay = ref_assay) %>% t()
  cell_states <- GetAssayData(slide, assay = state_assay) %>% t()
  
  list_mats <- map(names(cts) %>% set_names(), function(x) {
    
    mask_state_matrix(cell_props = cell_props, 
                      cell_states = cell_states,
                      cell_type = x,
                      cell_type_prop = cts[x])
    
  })
  
  mod_states <- reduce(list_mats, cbind)
  
  return(mod_states)
  
}

# Main

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
out_folder = "./results/state_scores/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples,
                    out_file = paste0(out_folder, 
                                      sample,"_state_map.pdf"))

# Create function to mask spots where the cell may not be present ---------------------------------------------
cts <- set_names(c(0.125, 0.125, 0.125),
                 c("cardiomyocyte", "endothelial", 
                   "fibroblast"))

walk2(visium_df$visium_file, visium_df$out_file, function(visium_file, out_file) { 
  
  print(visium_file)
  
  visium_slide <- readRDS(visium_file) 
  
  visium_slide[["cell_states"]] <- CreateAssayObject(data = t(filter_state_assay(cts = cts, 
                                                                                 slide = visium_slide)))
  DefaultAssay(visium_slide) <- "cell_states"
  
  pdf(file = out_file, height = 11, width = 14)
  
  walk(names(cts), function(cell_type) {
    
    features_ix <- grepl(cell_type, rownames(visium_slide))
    features <- rownames(visium_slide)[features_ix]
    
    states_plts <- map(features, function(spec_f) {
      
      SpatialFeaturePlot(visium_slide, 
                         features = spec_f) 
      
    })
    
    niche_plots <- plot_grid(plotlist= states_plts, ncol = 3)
    
    plot(niche_plots)
  })
  
  dev.off()

})


