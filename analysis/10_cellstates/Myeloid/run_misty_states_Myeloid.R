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
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# State descriptions 
Fib_state_list <- list("group_1" = c("Fib-state1", "Fib-state2", "Fib-state3"),
                   "group_2" = c("Fib-state2", "Fib-state5", "Fib-state6", "Fib-state7"),
                   "group_3" = c("Fib-state1", "Fib-state2", "Fib-state3", "Fib-state5"))

Myeloid_state_list <- list("group_1" = c("Myeloid-state0", "Myeloid-state2", "Myeloid-state3"),
                           "group_2" = c("Myeloid-state0", "Myeloid-state3", "Myeloid-state5", "Myeloid-state6"),
                           "group_3" = c("Myeloid-state0", "Myeloid-state3", "Myeloid-state2", "Myeloid-state4"))

# Cell2location scores - complete
assay_label <- "c2l"
ct <- "Myeloid"

run_state_ppline <- function(ct, assay_label){
  
  print(ct)
  
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
    all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states"))
    # Make the choice based on patient group
    pat_group <- sample_dict %>%
      dplyr::filter(sample_id == slide_id) %>%
      pull(patient_group)
    
    # Get Myeloid features
    useful_ct_states <- Myeloid_state_list[[pat_group]]
    ct_states <- all_ct_states[all_ct_states %in% useful_ct_states]
    
    # Get Fib features
    useful_ct_states_f <- Fib_state_list[[pat_group]]
    ct_states_f <- all_ct_states[all_ct_states %in% useful_ct_states_f]
    
    misty_out_alias <- paste0("./results/state_structure/", ct,"/", pat_group, "/mstate_")
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = "cell_states",
                        "intra_cts" = assay_label,
                        "para_states" = "cell_states",
                        "para_cts" = assay_label,
                        "intra_fib" = "cell_states",
                        "para_fib" = "cell_states")
    
    # Define features of each view ------------
    view_features <- list("main" = ct_states,
                          "intra_cts" = cts,
                          "para_states" = ct_states,
                          "para_cts" = cts,
                          "intra_fib" = ct_states_f,
                          "para_fib" = ct_states_f) #Use all
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "intra_cts" = "intra",
                       "para_states" = "para",
                       "para_cts" = "para",
                       "intra_fib" = "intra",
                       "para_fib" = "para") #Use all
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "intra_cts" = NULL,
                        "para_states" = 15, #Use all
                        "para_cts" = 15,
                        "intra_fib" = NULL,
                        "para_fib" = 15) 
    
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
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_cts_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_states_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_fib", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_fib_15", cutoff = 0)
    
    dev.off()
    
    return(misty_out)
    
  })
}

#Main

run_state_ppline(ct = ct,assay_label = assay_label)

#Now group them per folder

patient_annotation <- sample_dict

performance_all_misty = function(misty_out_folder, r2_filter = 0.1) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .))
  
  misty_res <- collect_results(misty_outs)
  
  model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                  grepl("multi.R2", measure)) %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .))
  
  importances_filtered <- misty_res$importances %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .))
  
  best_performers <- model_performance %>% 
    dplyr::filter(measure == "multi.R2") %>% 
    group_by(target) %>% 
    dplyr::filter(value > r2_filter) %>%
    arrange(target) 
  
  
  selected_importances <- best_performers %>%
    left_join(importances_filtered, 
              by = c("sample", 
                     "target" = "Target"))
  
  median_importances <- selected_importances %>%
    group_by(target, view, Predictor) %>%
    summarise(median_imp = median(Importance)) %>%
    ungroup() %>%
    group_by(view) %>%
    nest()
  
  pdf(file = paste0(misty_out_folder,"/summary_plot_mean.pdf"))
  
  walk2(median_importances$view, 
        median_importances$data, function(v, dat){
          
          imp_plot <- dat %>%
            ggplot(aes(x = target, y = Predictor, fill = median_imp)) +
            geom_tile() +
            scale_fill_gradient2(high = "#8DA0CB", 
                                 midpoint = 0,
                                 low = "white",
                                 na.value = "grey") +
            ggplot2::coord_equal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(v)
          
          print(imp_plot)
        })
  
  dev.off()
  
  best_performers %>%
    write_csv(file = paste0(misty_out_folder,"/best_performers.csv"))
  
}

performance_all_misty("./results/state_structure/Myeloid/group_1/")
performance_all_misty("./results/state_structure/Myeloid/group_2/")
performance_all_misty("./results/state_structure/Myeloid/group_3/")

####

summarize_all_misty = function(misty_out_folder, r2_filter = 0.2) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .))
  
  misty_res <- collect_results(misty_outs)
  
  model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                  grepl("multi.R2", measure)) %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .)) %>%
    left_join(patient_annotation, by = c("sample" = "sample_id"))
  
  mean_multi <- model_performance %>% 
    dplyr::filter(measure == "multi.R2") %>% 
    group_by(sample,patient_group) %>% 
    summarise(meanR2 = mean(value)) %>% 
    arrange(-meanR2)
  
  selected_samples <- dplyr::filter(mean_multi, meanR2 >= r2_filter) %>%
    pull(sample)
  
  selected_files <- paramdf %>%
    dplyr::filter(sample %in% selected_samples)
  
  misty_res <- collect_results(selected_files$misty_out_file)
  
  pdf(file = paste0(misty_out_folder,"/summary_plot_mean.pdf"))
  
  mistyR::plot_improvement_stats(misty_res)
  mistyR::plot_view_contributions(misty_res)
  
  mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
  mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)
  
  mistyR::plot_interaction_heatmap(misty_res, "intra_cts", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res, "para_cts_15", cutoff = 0)
  
  mistyR::plot_interaction_heatmap(misty_res, "para_states_15", cutoff = 0)
  
  
  dev.off()
  
}

summarize_all_misty("./results/state_structure/Fib/group_1/")
summarize_all_misty("./results/state_structure/Fib/group_2/")
summarize_all_misty("./results/state_structure/Fib/group_3/")

misty_out_folder = "./results/state_structure/Fib/group_1/"


