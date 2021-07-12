# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run MISTy for spatial interactions between slides
library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_pplne_utils.R")

# Pipeline definition:
run_colocalization <- function(slide, assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./visium_results_manuscript/colocalization/mjr_") {
  
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
                      "para" = 15)
  
  misty_out <- paste0(misty_out_alias, 
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
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Cell2location proportions - complete
assay_label <- "c2l_major_props"

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
                     assay = assay,
                     misty_out_alias =  "./visium_results_manuscript/colocalization/all_")
  
  return(mout)
  
})

misty_res <- mistyR::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_allc2lprops.pdf")

print(misty_res$improvements %>%
    dplyr::filter(grepl("R2", measure) &
                    !grepl("p", measure) &
                    !grepl("gain", measure)) %>%
  ggplot(aes(x = target, color = measure, y = value)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "juxta_5", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "para_15", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "para_15", cutoff = 0.5)

dev.off()

# Cell2location proportions - filtered
assay_label <- "c2l_major_props"

misty_outs <- map(slide_files, function(slide_file){
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  useful_features <- useful_features[!grepl("cardiomyocytes", useful_features)]
  
  mout <- run_colocalization(slide = slide,
                             useful_features = useful_features,
                             out_label = slide_id,
                             assay = assay,
                             misty_out_alias =  "./visium_results_manuscript/colocalization/wocardio_")
  
  return(mout)
  
})

misty_res <- mistyR::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_wocardio_c2lprops.pdf")

print(misty_res$improvements %>%
        dplyr::filter(grepl("R2", measure) &
                        !grepl("p", measure) &
                        !grepl("gain", measure)) %>%
        ggplot(aes(x = target, color = measure, y = value)) +
        geom_violin() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "juxta_5", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "para_15", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "para_15", cutoff = 0.5)

dev.off()

# Cell2location  - complete
assay_label <- "c2l_major"

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
                             assay = assay,
                             misty_out_alias =  "./visium_results_manuscript/colocalization/all_")
  
  return(mout)
  
})

misty_res <- mistyR::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_allc2l.pdf")

print(misty_res$improvements %>%
        dplyr::filter(grepl("R2", measure) &
                        !grepl("p", measure) &
                        !grepl("gain", measure)) %>%
        ggplot(aes(x = target, color = measure, y = value)) +
        geom_violin() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "juxta_5", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "para_15", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "para_15", cutoff = 0.5)

dev.off()

# Cell2location proportions - filtered
assay_label <- "c2l_major"

misty_outs <- map(slide_files, function(slide_file){
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  assay <- assay_label
  DefaultAssay(slide) <- assay
  useful_features <- rownames(slide)
  useful_features <- useful_features[!grepl("cardiomyocytes", useful_features)]
  
  mout <- run_colocalization(slide = slide,
                             useful_features = useful_features,
                             out_label = slide_id,
                             assay = assay,
                             misty_out_alias =  "./visium_results_manuscript/colocalization/wocardio_")
  
  return(mout)
  
})

misty_res <- mistyR::collect_results(unlist(misty_outs))

pdf("./visium_results_manuscript/colocalization/misty_colocalization_wocardio_c2l.pdf")

print(misty_res$improvements %>%
        dplyr::filter(grepl("R2", measure) &
                        !grepl("p", measure) &
                        !grepl("gain", measure)) %>%
        ggplot(aes(x = target, color = measure, y = value)) +
        geom_violin() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))

mistyR::plot_improvement_stats(misty_res)
mistyR::plot_view_contributions(misty_res)

mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "juxta_5", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "juxta_5", cutoff = 0.5)

mistyR::plot_interaction_heatmap(misty_res, "para_15", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "para_15", cutoff = 0.5)

dev.off()
