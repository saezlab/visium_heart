# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Estimate colocalization with ct scores in
#' individual slides.
#' 
#' Fits NMF for each data set
#' 
library(Seurat)
library(tidyverse)
library(NNLM)

# Define a single slide analysis

colocalize_cts <- function(ct_mat, 
                           slide_file, 
                           out_file,
                           alias = "spotlight",
                           alias_nmf = "spotlight_nmf",
                           k = 5) {
  print("Analyzing...")
  print(slide_file)
  
  deconv_mat <- readRDS(ct_mat)
  slide <- readRDS(slide_file)
  slide[[alias]] <- CreateAssayObject(data = t(deconv_mat))
  
  decomp <- nnmf(deconv_mat, 
                 k = k, 
                 max.iter = 500)
  
  # Making objects compatible with Seurat structure
  factor_weights <- decomp$W
  colnames(factor_weights) <- paste0("factor_",1:ncol(factor_weights))
  
  factor_loadings <- decomp$H
  rownames(factor_loadings) <- paste0("factor_",1:nrow(factor_loadings))
  
  slide[[paste0("red_",alias_nmf)]] <- CreateDimReducObject(embeddings = factor_weights, 
                                             key = "NMF_", 
                                             assay = alias,
                                             loadings = t(factor_loadings))
  
  slide[[alias_nmf]] <- CreateAssayObject(data = t(factor_weights))
  
  saveRDS(slide, file = out_file)
  
  # Plotting ct_scores ------------------------------------------------
  #DefaultAssay(slide) <- alias_nmf
  #ctplts <- SpatialFeaturePlot(slide, features = rownames(slide))
  
  # Plotting factors ------------------------------------------------
  DefaultAssay(slide) <- alias_nmf
  factorplts <- SpatialFeaturePlot(slide, features = rownames(slide))
  
  # Plotting loadings ------------------------------------------------
  factor_loadings <- factor_loadings %>%
    as.data.frame() %>%
    rownames_to_column("factor") %>%
    pivot_longer(-factor, 
                 names_to = "feature", 
                 values_to = "loading")
  
  loadingplts <- ggplot(factor_loadings, aes(fill = loading, x = factor, y = feature)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list("factors_plts" = factorplts, 
              "loading_plts" = loadingplts))
  
}

# MAIN -------------------------
# Get slides and deconvolution matrices from SPOTlight
deconv_mats_folder <- "./visium_results_manuscript/deconvolution/spotlight_mjr/"
slide_files_folder <- "./visium_results_manuscript/processed_visium/"
out_slides_folder <- "./visium_results_manuscript/processed_visium_revisions/"
deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
slide_files <- list.files(slide_files_folder)
data_df <- tibble(ct_mat = paste0(deconv_mats_folder, sort(deconv_mats)),
                  slide_file = paste0(slide_files_folder, sort(slide_files)),
                  out_file = paste0(out_slides_folder, sort(slide_files)))


walk(c(4,5,6), function(k) { 
  
  print(k)
  
  all_plts <- pmap(data_df, colocalize_cts, k = k)
  
  pdf(file = paste0("./visium_results_manuscript/deconvolution/spotlight_qc_mjr/single_colocalization_",k,".pdf"), 
      height = 12, 
      width = 12)
  
  walk(all_plts, function(x) {
    
    walk(x, print)
    
  })
  
  dev.off()
  
})

# Get slides and deconvolution matrices from cell2location
deconv_mats_folder <- "./visium_results_manuscript/deconvolution/c2l/location_models/density_tables_rds/"
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
out_slides_folder <- "./visium_results_manuscript/processed_visium_revisions/"
deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
slide_files <- list.files(slide_files_folder)
data_df <- tibble(ct_mat = paste0(deconv_mats_folder, sort(deconv_mats)),
                  slide_file = paste0(slide_files_folder, sort(slide_files)),
                  out_file = paste0(out_slides_folder, sort(slide_files)))


walk(c(3, 4, 5, 6), function(k) { 
  
  print(k)
  
  all_plts <- pmap(data_df, 
                   colocalize_cts, 
                   k = k, 
                   alias = "c2l",
                   alias_nmf = "c2l_nmf")
  
  pdf(file = paste0("./visium_results_manuscript/deconvolution/c2l/colocalization/single_colocalization_",k,".pdf"), 
      height = 12, 
      width = 12)
  
  walk(all_plts, function(x) {
    
    walk(x, print)
    
  })
  
  dev.off()
  
})

# Get slides and deconvolution matrices from cell2location (states)
deconv_mats_folder <- "./visium_results_manuscript/deconvolution/c2l_states/location_models/density_tables_rds/"
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
out_slides_folder <- "./visium_results_manuscript/processed_visium_revisions/"
deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
slide_files <- list.files(slide_files_folder)
data_df <- tibble(ct_mat = paste0(deconv_mats_folder, sort(deconv_mats)),
                  slide_file = paste0(slide_files_folder, sort(slide_files)),
                  out_file = paste0(out_slides_folder, sort(slide_files)))


walk(c(4, 5, 6, 7), function(k) { 
  
  print(k)
  
  all_plts <- pmap(data_df, 
                   colocalize_cts, 
                   k = k, 
                   alias = "c2l_states",
                   alias_nmf = "c2l_nmf_states")
  
  pdf(file = paste0("./visium_results_manuscript/deconvolution/c2l_states/colocalization/single_colocalization_",k,".pdf"), 
      height = 12, 
      width = 12)
  
  walk(all_plts, function(x) {
    
    walk(x, print)
    
  })
  
  dev.off()
  
})

# Get slides and deconvolution matrices from cell2location (states)
deconv_mats_folder <- "./visium_results_manuscript/deconvolution/c2l_statescollapsed//location_models/density_tables_rds/"
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
out_slides_folder <- "./visium_results_manuscript/processed_visium_revisions/"
deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
slide_files <- list.files(slide_files_folder)
data_df <- tibble(ct_mat = paste0(deconv_mats_folder, sort(deconv_mats)),
                  slide_file = paste0(slide_files_folder, sort(slide_files)),
                  out_file = paste0(out_slides_folder, sort(slide_files)))


walk(c(4, 5, 6, 7), function(k) { 
  
  print(k)
  
  all_plts <- pmap(data_df, 
                   colocalize_cts, 
                   k = k, 
                   alias = "c2l_statesred",
                   alias_nmf = "c2l_nmf_statesred")
  
  pdf(file = paste0("./visium_results_manuscript/deconvolution/c2l_statescollapsed/colocalization/single_colocalization_",k,".pdf"), 
      height = 12, 
      width = 12)
  
  walk(all_plts, function(x) {
    
    walk(x, print)
    
  })
  
  dev.off()
  
})


