# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Test a simplified version of the spatial analysis of interacting cells of interest

library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")

future::plan(future::multisession)

#This generates an or statement of cts
get_ct_spots <- function(slide, ct, filter_prop = 0.05) {
  
  all_ids <- map(ct, function(x) {
    c_props <- GetAssayData(slide, assay = "c2l_props")[x,]
    names(c_props)[c_props >= filter_prop]
  })
  
  purrr::reduce(all_ids, intersect)
  
}

# Main ------------------------------------------------------------------------
# Getting sample annotations --------------------------------------------------
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Targets
target_list <- c("Fib-state0","Fib-state1", 
                 "Fib-state2", "Fib-state3", 
                 "Fib-state4")

predictors <- c("vSMCs-state0", "vSMCs-state1")

# Cell2location scores - complete
assay_label <- "c2l"
ct <- c("Fib","vSMCs")

# Main pipeline
run_state_ppline <- function(ct, assay_label){
  
  misty_outs <- map(slide_files, function(slide_file){
    
    print(slide_file)
    
    # Read spatial transcriptomics data
    slide_id <- gsub("[.]rds", "", slide_file)
    slide <- readRDS(paste0(slide_files_folder, slide_file))
    
    # Get useful spots
    useful_spots <- get_ct_spots(slide, ct)
    
    # Get ct states that are useful
    all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states"))
    predictor_states <- predictors[predictors %in% all_ct_states]
    
    misty_out_alias <- paste0("./results/state_structure/", 
                              paste0(ct,collapse = "_"), 
                              "/mstate_")
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = "cell_states",
                        "para_states" = "cell_states",
                        "intra_others" = "cell_states",
                        "para_others" = "cell_states")
    
    # Define features of each view ------------
    view_features <- list("main" = target_list,
                          "para_states" = target_list,
                          "intra_others" = predictor_states,
                          "para_others" = predictor_states) #Use all
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "para_states" = "para",
                       "intra_others" = "intra",
                       "para_others" = "para") #Use all
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "para_states" = 5, #Use all
                        "intra_others" = NULL,
                        "para_others" = 5) 
    
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
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_states_5", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "intra_others", cutoff = 0)
    
    mistyR::plot_interaction_heatmap(misty_res_slide, "para_others_5", cutoff = 0)
    
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
             gsub("_c2l","", .)) %>%
    left_join(sample_dict, by = c("sample" = "sample_id"))
  
  cors_filtered <- paramdf %>%
    dplyr::select(sample, cor_res) %>%
    unnest() %>%
    dplyr::rename("correlation" = value) %>%
    dplyr::mutate(feature_a = gsub("-", ".", feature_a),
                  feature_b = gsub("-", ".", feature_b)) %>%
    left_join(sample_dict, by = c("sample" = "sample_id"))
  
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
  
  
  mean_importances <- selected_importances %>%
    group_by(target, view, Predictor, patient_group) %>%
    summarise(mean_imp = mean(Importance)) %>%
    ungroup() %>%
    group_by(view) %>%
    nest()
  
  mean_cors <- selected_cors %>%
    group_by(target, feature_b, patient_group) %>%
    summarise(mean_cor = mean(correlation)) %>%
    ungroup() %>%
    dplyr::filter(target != feature_b)
  
  pdf(file = paste0(misty_out_folder,"/summary_plot_mean.pdf"))
  
  walk2(mean_importances$view, 
        mean_importances$data, function(v, dat){
          
          imp_plot <- dat %>%
            ggplot(aes(x = target, y = Predictor, fill = mean_imp)) +
            geom_tile() +
            scale_fill_gradient2(high = "#8DA0CB", 
                                 midpoint = 0,
                                 low = "white",
                                 na.value = "grey") +
            ggplot2::coord_equal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(v) +
            facet_wrap(.~patient_group)
          
          print(imp_plot)
        })
  
  ct_cors_plt <- ggplot(mean_cors, 
                        aes(x = target, y = feature_b, fill = mean_cor)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("Mean correlation of \n states in best predictors") +
    facet_wrap(.~patient_group)
  
  plot(ct_cors_plt)
  
  dev.off()
  
  best_performers %>%
    write_csv(file = paste0(misty_out_folder,"/best_performers.csv"))
  
}

# Main
run_state_ppline(ct = ct,assay_label = assay_label)
misty_out_folder <- "./results/state_structure/Fib_vSMCs"
performance_all_misty(misty_out_folder, r2_filter = 0.1)

best_performers <- read_csv("./results/state_structure/Fib_vSMCs/best_performers.csv") %>%
  left_join(sample_dict, by = c("sample" = "sample_id"))

best_performers %>%
  ggplot(aes(x = patient_group, y = value)) +
  geom_boxplot() +
  facet_wrap(.~target)
