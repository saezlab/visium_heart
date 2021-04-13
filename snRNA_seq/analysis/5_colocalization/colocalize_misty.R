# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run MISTy for spatial interactions between slides
library(tidyverse)
library(Seurat)
library(MISTy)
source("./visium/visiumtools/pipelines.R")

# Pipeline definition:
run_colocalization <- function(slide, assay, useful_features, out_label) {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 5,
                      "para" = 10)
  
  misty_out <- paste0("./visium_results_manuscript/colocalization/mjr_", 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}

# Main -----------------------------------
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Spotlight
assay_label <- "spotlight"

misty_outs <- map(slide_files, function(slide_file){
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  useful_features <- useful_features[!grepl("state-cardiomyocytes", useful_features)]
  
  mout <- run_colocalization(slide = slide,
                     useful_features = useful_features,
                     out_label = slide_id,
                     assay = assay)
  
  return(mout)
  
})

misty_res <- MISTy::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_spotlight.pdf")
MISTy::plot_improvement_stats(misty_res)
MISTy::plot_view_contributions(misty_res)
MISTy::plot_interaction_heatmap(misty_res, "intra", cutoff = 0.5)
MISTy::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0.5)
MISTy::plot_interaction_heatmap(misty_res, "para_10", cutoff = 0.5)
dev.off()

# cell2location

assay_label <- "c2l"

misty_outs <- map(slide_files, function(slide_file){
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  mout <- run_colocalization(slide = slide,
                             useful_features = useful_features,
                             out_label = slide_id,
                             assay = assay)
  
  return(mout)
  
})

misty_res <- MISTy::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_c2l.pdf")
MISTy::plot_improvement_stats(misty_res)
MISTy::plot_view_contributions(misty_res)
MISTy::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
MISTy::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
MISTy::plot_interaction_heatmap(misty_res, "para_10", cutoff = 0)
dev.off()

