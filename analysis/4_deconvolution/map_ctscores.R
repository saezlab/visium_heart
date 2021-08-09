library(Seurat)
library(tidyverse)

deconv_mats_path <- "./visium_results_manuscript/deconvolution/spotlight_mjr/"
slide_files_path <- "./visium_results_manuscript/processed_visium/"

deconv_mats <- list.files(deconv_mats_path)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]

slide_files <- list.files(slide_files_path)

data_df <- tibble(ct_mat = sort(deconv_mats),
                  slide = sort(slide_files))

data_df$ct_mat = paste0(deconv_mats_path, 
                        data_df$ct_mat)

data_df$slide = paste0(slide_files_path, 
                        data_df$slide)

printCT <- function(ct_mat, slide_file, features = c("state-CM1", 
                                                     "state-CM2",
                                                     "state-CM3",
                                                     "state-Damaged-CM")) {
  
  deconv_mat <- readRDS(ct_mat)
  slide <- readRDS(slide_file)
  slide[['spotlight']] <- CreateAssayObject(data = t(deconv_mat))
  DefaultAssay(slide) <- "spotlight"
  cell_states_plots <- SpatialFeaturePlot(slide, features, ncol = 3)
  
}

#CMs

data_df <- data_df %>%
  dplyr::mutate(cm_plts = map2(ct_mat,slide, printCT))

pdf("./visium_results_manuscript/deconvolution/spotlight_cms.pdf", height = 8, width = 10)

walk(data_df$cm_plts, function(x) { 
  
  print(x)
  
  })

dev.off()


#FBs

data_df <- data_df %>%
  dplyr::mutate(cm_plts = map2(ct_mat,slide, printCT, features = c("state-Fib1",
                                                                   "state-Fib2",
                                                                   "state-Fib3SCARA5",
                                                                   "state-Fib4",
                                                                   "state-Myofib")))

pdf("./visium_results_manuscript/deconvolution/spotlight_fbs.pdf", height = 8, width = 10)

walk(data_df$cm_plts, function(x) { 
  
  print(x)
  
})

dev.off()

#Rest

data_df <- data_df %>%
  dplyr::mutate(cm_plts = map2(ct_mat,slide, printCT, features = c("state-adipocytes",
                                                                   "state-cardiomyocytes",
                                                                   "state-fibroblasts",
                                                                   "state-endothelial-cells",
                                                                   "state-lymphatic-endo",
                                                                   "state-macrophages",
                                                                   "state-pericytes",
                                                                   "state-t-cells")))

pdf("./visium_results_manuscript/deconvolution/spotlight_mjr.pdf", height = 8, width = 10)

walk(data_df$cm_plts, function(x) { 
  
  print(x)
  
})

dev.off()


