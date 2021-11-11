# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Test a simplified version of the spatial analysis of interacting cells of interest

library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")

future::plan(future::multisession)

#This generates an or statement of cts
get_ct_spots <- function(slide, ct) {
  rownames(slide@meta.data)[slide@meta.data[,paste0(ct,"_flag")] == 1]
}

# Main ------------------------------------------------------------------------

# Getting sample annotations --------------------------------------------------
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Targets
target_list <- c("Fib-state0","Fib-state1", "Fib-state2", "Fib-state3", "Fib-state4")
predictors <- c("Myeloid-state0", "Myeloid-state1", "Myeloid-state2", 
                "Myeloid-state3", "vSMCs-state0", "vSMCs-state1")

# Cell2location scores - complete
assay_label <- "c2l"
ct <- "Fib"

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
    predictor_states <- predictors[predictors %in% all_ct_states]
    
    misty_out_alias <- paste0("./results/state_structure/", ct, "/mstate_")
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = "cell_states",
                        "intra_cts" = assay_label,
                        "para_states" = "cell_states",
                        "para_cts" = assay_label,
                        "intra_others" = "cell_states",
                        "para_others" = "cell_states")
    
    # Define features of each view ------------
    view_features <- list("main" = target_list,
                          "intra_cts" = cts,
                          "para_states" = target_list,
                          "para_cts" = cts,
                          "intra_others" = predictor_states,
                          "para_others" = predictor_states) #Use all
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "intra_cts" = "intra",
                       "para_states" = "para",
                       "para_cts" = "para",
                       "intra_others" = "intra",
                       "para_others" = "para") #Use all
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "intra_cts" = NULL,
                        "para_states" = 15, #Use all
                        "para_cts" = 15,
                        "intra_others" = NULL,
                        "para_others" = 15) 
    
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
    
    # Plot correlations for the sake of completion
    state_data <- GetAssayData(slide, assay = "cell_states")[c(predictor_states, target_list),]
    
    cor_data_states <- cor(t(state_data)) %>%
      as.data.frame() %>%
      rownames_to_column("feature_a") %>%
      pivot_longer(-feature_a, names_to = "feature_b")
    
    write_csv(cor_data_states, file = paste0(plot_folder,"/state_cor.csv"))
    
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
    
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
    mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_cts", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_cts_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_states_15", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_others", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_others_15", cutoff = 0)
    
    cor_data_states_plt <- cor_data_states %>%
      ggplot(aes(x = feature_a, y = feature_b, fill = value)) +
      geom_tile() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_gradient2(midpoint = 0)
    
    plot(cor_data_states_plt)
    
    dev.off()
    
    return(misty_out)
    
  })
}

# This is what groups and selects the relevant slides

performance_all_misty = function(misty_out_folder, r2_filter = 0.1) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]]) %>%
             gsub("_c2l","", .)) %>%
    mutate(cor_res = paste0(misty_out_file, "/plots/state_cor.csv")) %>%
    mutate(cor_res = map(cor_res, read_csv))
  
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
  
  cors_filtered <- paramdf %>%
    dplyr::select(sample, cor_res) %>%
    unnest() %>%
    dplyr::rename("correlation" = value) %>%
    dplyr::mutate(feature_a = gsub("-", ".", feature_a),
                  feature_b = gsub("-", ".", feature_b))
  
  best_performers <- model_performance %>% 
    dplyr::filter(measure == "multi.R2") %>% 
    group_by(target) %>% 
    dplyr::filter(value > r2_filter) %>%
    arrange(target) 
  
  selected_importances <- best_performers %>%
    left_join(importances_filtered, 
              by = c("sample", 
                     "target" = "Target"))
  
  selected_cors <- best_performers %>%
    left_join(cors_filtered, 
              by = c("sample", 
                     "target" = "feature_a"))
     
  
  median_importances <- selected_importances %>%
    group_by(target, view, Predictor) %>%
    summarise(median_imp = median(Importance)) %>%
    ungroup() %>%
    group_by(view) %>%
    nest()
  
  median_cors <- selected_cors %>%
    group_by(target, feature_b) %>%
    summarise(median_cor = median(correlation)) %>%
    ungroup()
  
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
  
  ggplot(aes(x = target, y = feature_b, fill = median_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Median correlation of \n states in best predictors")
  
  dev.off()
  
  best_performers %>%
    write_csv(file = paste0(misty_out_folder,"/best_performers.csv"))
  
}



run_state_ppline(ct = ct,assay_label = assay_label)
misty_out_folder <- "./results/state_structure/Fib"






