library(Seurat)
library(tidyverse)

deconv_mats <- list.files("./visium_results_manuscript/deconvolution/spotlight/")
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]

slide_files <- list.files("./visium_results_manuscript/processed_visium/")

data_df <- tibble(ct_mat = sort(deconv_mats),
                  slide = sort(slide_files))

data_df$ct_mat = paste0("./visium_results_manuscript/deconvolution/spotlight/", 
                        data_df$ct_mat)

data_df$slide = paste0("./visium_results_manuscript/processed_visium/", 
                        data_df$slide)

printCT <- function(ct_mat, slide_file, features = c("state-CM1", 
                                                     "state-CM2",
                                                     "state-CM3",
                                                     "state-Damaged-CM")) {
  
  deconv_mat <- readRDS(ct_mat)
  slide <- readRDS(slide_file)
  slide[['spotlight']] <- CreateAssayObject(data = t(deconv_mat))
  DefaultAssay(slide) <- "spotlight"
  cell_states_plots <- SpatialFeaturePlot(slide, features, ncol = 2)
  
}


data_df <- data_df %>%
  dplyr::mutate(cm_plts = map2(ct_mat,slide, printCT))





