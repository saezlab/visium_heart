# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we will run MISTy to explain the cell-state scores of different cell-types
#' First, we require to define meta variables per spot that flag spots with
#' enough cell-types
#' 
#' Once the spots, MISTy will be run using this structure:
#' 
#' 1) Intra Main - states
#' 2) 
#' 2) Para - cell2location cell scores (without cell of states)
#' 
library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")

future::plan(future::multisession)

get_ct_spots <- function(slide, ct) {
  rownames(slide@meta.data)[slide@meta.data[,paste0(ct,"_flag")] == 1]
}

# Main ------------------------------------------------------------------------

# Getting sample annotations --------------------------------------------------

sample_dict <- read.table("./markers/visium_annotations_ext.txt", 
                          sep = "\t", header = T)

slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Cell2location scores - complete
assay_label <- "c2l"
ct <- "Fib"

run_state_ppline <- function(ct, assay_label){
  
  print(ct)
  
  misty_out_alias <- paste0("./results/state_structure/", ct, "/mstate_")
  
  misty_outs <- map(slide_files, function(slide_file){
    
    print(slide_file)
    
    # Read spatial transcriptomics data
    slide_id <- gsub("[.]rds", "", slide_file)
    slide <- readRDS(paste0(slide_files_folder, slide_file))
    
    # Get useful spots
    useful_spots <- get_ct_spots(slide, ct)
    
    # Get ct features that are not the cell of interest
    cts <- rownames(GetAssayData(slide, assay = assay_label))
    cts <- cts[!grepl(ct, cts)]
    
    # Get ct states that are useful
    ct_states <- rownames(GetAssayData(slide, assay = "cell_states"))
    ct_states <- ct_states[grepl(ct, ct_states)]
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = "cell_states",
                        "intra_cts" = assay_label,
                        "intra_progeny" = "progeny",
                        "para_states" = "cell_states",
                        "para_cts" = assay_label,
                        "para_progeny" = "progeny")
    
    # Define features of each view ------------
    view_features <- list("main" = ct_states,
                          "intra_cts" = cts,
                          "intra_progeny" = NULL, #Use all,
                          "para_states" = ct_states,
                          "para_cts" = cts,
                          "para_progeny" = NULL) #Use all
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "intra_cts" = "intra",
                       "intra_progeny" = "intra", #Use all
                       "para_states" = "para",
                       "para_cts" = "para",
                       "para_progeny" = "para") #Use all
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "intra_cts" = NULL,
                        "intra_progeny" = NULL, #Use all
                        "para_states" = 15, #Use all
                        "para_cts" = 15,
                        "para_progeny" = 15) 
    
    misty_out <- paste0(misty_out_alias, 
                        slide_id, "_", assay_label)
    
    run_misty_seurat(visium.slide = slide,
                     view.assays = view_assays,
                     view.features = view_features,
                     view.types = view_types,
                     view.params = view_params,
                     spot.ids = useful_spots,
                     out.alias = misty_out)
    
    misty_res_slide <- collect_results(misty_out)
    
    plot_folder <- paste0(misty_out, "/plots")
    
    system(paste0("mkdir ", plot_folder))
    
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
    
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
    mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_cts", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_progeny", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_cts_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_states_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_progeny_15", cutoff = 0)
    
    dev.off()
    
    return(misty_out)
    
  })
}

# 

run_state_ppline(ct = "Fib",
                 assay_label = assay_label)

paramdf <- tibble(cts = c("fibroblast", "cardiomyocyte", "endothelial"),
                  assay_labels = "c2l")

walk2(paramdf$cts, paramdf$assay_labels, run_state_ppline)

