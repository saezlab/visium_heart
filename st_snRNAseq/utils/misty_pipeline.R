# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Functions used in pipelines
library(tidyverse)
library(Seurat)
library(mistyR)

#This generates a region of interest based on the proportion
#of a collection of cell-types
get_ct_spots <- function(slide, ct, filter_prop = 0.1) {
  
  all_ids <- map(ct, function(x) {
    c_props <- GetAssayData(slide, assay = "c2l_props")[x,]
    names(c_props)[c_props >= filter_prop]
  })
  
  purrr::reduce(all_ids, intersect)
  
}

#This generates a mask of negative scores
positive_states <- function(slide, assay = "cell_states") {
  
  state_mats <- GetAssayData(slide, assay = assay)
  
  state_mats[state_mats < 0] <- 0
  
  slide[["cell_states_pos"]] <- CreateAssayObject(data = state_mats)
  
  return(slide)
}

#This generates a mask for non-informative spots based on compositions
# by_prop: TRUE
filter_states <- function(slide, 
                          by_prop = T, 
                          prop_thrsh = 0.1,
                          assay = "cell_states_pos") {
  
  # Get the cell-types associated with the states
  # As a name they have ct-state_id
  state_cms <- rownames(GetAssayData(slide, assay = assay)) %>%
    strsplit(., "-") %>%
    map_chr(., ~.x[[1]]) %>%
    unique() %>%
    set_names()
  
  state_mats <- map(state_cms, function(ct) {
    print(ct)
    
    prop_vector <- GetAssayData(slide, assay = "c2l_props")[ct,]
    
    #If you don't want to consider pure proportions use:
    if(!by_prop){
      prop_vector <- ifelse(prop_vector >= prop_thrsh, 1, 0)
    }
    
    state_mat <- GetAssayData(slide, assay = assay)
    
    ix <-  rownames(state_mat)[grepl(ct, rownames(state_mat))]
    state_mat <- state_mat[ix, ]
    
    apply(state_mat, 1, function(x) {
      x * prop_vector
    }) %>%
      t()
    
  })
  
  state_mats <- purrr::reduce(state_mats, rbind)
  
  slide[["cell_states_pos"]] <- CreateAssayObject(data = state_mats)
  
  return(slide)
}


# Main pipeline
# This models interactions between cell-states in ROIs

run_state_ppline <- function(ROI_ct, #cts defining the region of interest
                             ROI_prop = 0.1, # cut-off that defines ROI (intersection of all spots)
                             mask_by_prop = F, # this is a mask based on proportion of cells
                             mask_threshold = 0.1, # This is a simpler mask based on ROI 
                             folder_label, # how to call the MISTy res folder
                             targets,
                             predictors,
                             target_assay = "cell_states_pos",
                             join_trgts_pred = T, # If target and predictors are modeled together, then we only do cell_states_pos
                             state_origin = "cell_states"){
  
  misty_outs <- map(slide_files, function(slide_file){
    
    print(slide_file)
    
    # Read spatial transcriptomics data and transform states to be useful for modeling
    slide_id <- gsub("[.]rds", "", slide_file)
    slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
      positive_states(., assay = state_origin)
    
    # This modifications stay in cell_states_pos
    if(mask_by_prop) {
      
      slide <- filter_states(slide = slide,
                             by_prop = T,
                             prop_thrsh = NULL)
    } else {
      slide <- filter_states(slide = slide,
                             by_prop = F,
                             prop_thrsh = mask_threshold)
    }
    
    # Get useful spots - ROI
    useful_spots <- get_ct_spots(slide, ct = ROI_ct, filter_prop = ROI_prop)
    
    if(length(useful_spots) < 25) {
      
      return(NULL)
      
    } 
    
    # Get ct states that are useful
    
    # Here predictors and predicted values are considered together 
    # (IMPORTANCES of predictors may drop)
    if(join_trgts_pred) {
      
      model_states <- c(targets, predictors)
      all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states_pos"))
      all_ct_states <- all_ct_states[all_ct_states %in% model_states]
      
      # Define pipeline
      
      # Define assay of each view ---------------
      view_assays <- list("main" = "cell_states_pos",
                          "para_states" = "cell_states_pos")
      
      # Define features of each view ------------
      view_features <- list("main" = all_ct_states,
                            "para_states" = all_ct_states) #Use all
      
      # Define spatial context of each view -----
      view_types <- list("main" = "intra",
                         "para_states" = "para") #Use all
      
      # Define additional parameters (l in case of paraview,
      # n of neighbors in case of juxta) --------
      view_params <- list("main" = NULL,
                          "para_states" = 5) 
      
    } else {
      # Here predictors and predicted values are considered separately
      # Easier to see the interaction between targets/predictors
      
      all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states_pos"))
      model_states <- c(targets, predictors)
      all_ct_states <- all_ct_states[all_ct_states %in% model_states]
      
      useful_targets <- all_ct_states[all_ct_states %in% targets]
      useful_predictors <- all_ct_states[all_ct_states %in% predictors]
      
      # Define pipeline
      
      # Define assay of each view ---------------
      view_assays <- list("main" = target_assay,
                          "para_states" = target_assay,
                          "intra_pred" = "cell_states_pos",
                          "para_pred" = "cell_states_pos")
      
      # Define features of each view ------------
      view_features <- list("main" = useful_targets,
                            "para_states" = useful_targets,
                            "intra_pred" = useful_predictors,
                            "para_pred" = useful_predictors)
      
      # Define spatial context of each view -----
      view_types <- list("main" = "intra",
                         "para_states" = "para",
                         "intra_pred" = "intra",
                         "para_pred" = "para")
      
      # Define additional parameters (l in case of paraview,
      # n of neighbors in case of juxta) --------
      view_params <- list("main" = NULL,
                          "para_states" = 5,
                          "intra_pred" = NULL,
                          "para_pred" = 5)
      
    }
    
    misty_out_alias <- paste0("./results/state_structure/", 
                              folder_label, 
                              "/mstate_")
    
    misty_out <- paste0(misty_out_alias, 
                        slide_id)
    
    run_misty_seurat(visium.slide = slide,
                     view.assays = view_assays,
                     view.features = view_features,
                     view.types = view_types,
                     view.params = view_params,
                     spot.ids = useful_spots,
                     out.alias = misty_out)
    
    misty_res_slide <- collect_results(misty_out)
    
    view_names <- misty_res_slide$importances$view %>% unique()
    
    plot_folder <- paste0(misty_out, "/plots")
    
    system(paste0("mkdir ", plot_folder))
    
    # Plot correlations for the sake of completion
    if(join_trgts_pred) {
      state_data <- GetAssayData(slide, assay = "cell_states_pos")[unique(all_ct_states), useful_spots]
    } else {
      state_data_A <- GetAssayData(slide, assay = target_assay)[useful_targets, useful_spots]
      state_data_B <- GetAssayData(slide, assay = "cell_states_pos")[useful_predictors, useful_spots]
      state_data <- rbind(state_data_A, state_data_B)
    }
    
    cor_data_states <- cor(t(state_data)) %>%
      as.data.frame() %>%
      rownames_to_column("feature_a") %>%
      pivot_longer(-feature_a, names_to = "feature_b")
    
    write_csv(cor_data_states, file = paste0(plot_folder,"/state_cor.csv"))
    
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
    
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    walk(view_names, function(v) {
      mistyR::plot_interaction_heatmap(misty_res_slide, v, cutoff = 0)
    })
    
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


# Main pipeline
run_state_ppline_ct <- function(ROI_ct,
                                ROI_prop = 0.1,
                                mask_by_prop = F,
                                mask_threshold = 0.1,
                                folder_label, 
                                targets,
                                target_assay = "cell_states_pos",
                                state_origin = "cell_states"){
  
  misty_outs <- map(slide_files, function(slide_file){
    
    print(slide_file)
    
    # Read spatial transcriptomics data and transform states to be useful for modeling
    slide_id <- gsub("[.]rds", "", slide_file)
    slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
      positive_states(.,assay = state_origin)
    
    if(mask_by_prop) {
      
      slide <- filter_states(slide = slide,
                             by_prop = T,
                             prop_thrsh = NULL)
    } else {
      slide <- filter_states(slide = slide,
                             by_prop = F,
                             prop_thrsh = mask_threshold)
    }
    
    # Get useful spots - ROI
    useful_spots <- get_ct_spots(slide, ct = ROI_ct, filter_prop = ROI_prop)
    
    if(length(useful_spots) < 25) {
      
      return(NULL)
      
    } 
    
    # Get predictor cell-types
    all_cts <- rownames(GetAssayData(slide, assay = "c2l"))
    all_cts <- all_cts[!(all_cts %in% c(ROI_ct, "prolif"))]
    
    # Get ct states that are useful
    
    # Here predictors and predicted values are considered separately
    # Easier to see the interaction between targets/predictors
    
    all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states_pos"))
    
    useful_targets <- all_ct_states[all_ct_states %in% targets]
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = target_assay,
                        "para_states" = "cell_states_pos",
                        "intra_pred" = "c2l",
                        "para_pred" = "c2l")
    
    # Define features of each view ------------
    view_features <- list("main" = useful_targets,
                          "para_states" = useful_targets,
                          "intra_pred" = all_cts,
                          "para_pred" = all_cts)
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "para_states" = "para",
                       "intra_pred" = "intra",
                       "para_pred" = "para")
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "para_states" = 5,
                        "intra_pred" = NULL,
                        "para_pred" = 5)
    
    misty_out_alias <- paste0("./results/state_structure/", 
                              folder_label, 
                              "/mstate_")
    
    misty_out <- paste0(misty_out_alias, 
                        slide_id)
    
    run_misty_seurat(visium.slide = slide,
                     view.assays = view_assays,
                     view.features = view_features,
                     view.types = view_types,
                     view.params = view_params,
                     spot.ids = useful_spots,
                     out.alias = misty_out)
    
    misty_res_slide <- collect_results(misty_out)
    
    view_names <- misty_res_slide$importances$view %>% unique()
    
    plot_folder <- paste0(misty_out, "/plots")
    
    system(paste0("mkdir ", plot_folder))
    
    # Plot correlations for the sake of completion
    state_data_A <- GetAssayData(slide, assay = target_assay)[useful_targets, useful_spots]
    state_data_B <- GetAssayData(slide, assay = "c2l")[all_cts, useful_spots]
    state_data <- rbind(state_data_A, state_data_B)
    
    cor_data_states <- cor(t(state_data)) %>%
      as.data.frame() %>%
      rownames_to_column("feature_a") %>%
      pivot_longer(-feature_a, names_to = "feature_b")
    
    write_csv(cor_data_states, file = paste0(plot_folder,"/state_cor.csv"))
    
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
    
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    walk(view_names, function(v) {
      mistyR::plot_interaction_heatmap(misty_res_slide, v, cutoff = 0)
    })
    
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


# Running pipelines with single cell-type
# 
#
run_state_ppline_ct_simple <- function(ROI_ct,
                                ROI_prop = 0.1,
                                mask_by_prop = F,
                                mask_threshold = 0.1,
                                folder_label, 
                                targets,
                                target_assay = "cell_states_pos",
                                state_origin = "cell_states"){
  
  misty_outs <- map(slide_files, function(slide_file){
    
    print(slide_file)
    
    # Read spatial transcriptomics data and transform states to be useful for modeling
    slide_id <- gsub("[.]rds", "", slide_file)
    slide <- readRDS(paste0(slide_files_folder, slide_file)) %>%
      positive_states(.,assay = state_origin)
    
    if(mask_by_prop) {
      
      slide <- filter_states(slide = slide,
                             by_prop = T,
                             prop_thrsh = NULL)
    } else {
      slide <- filter_states(slide = slide,
                             by_prop = F,
                             prop_thrsh = mask_threshold)
    }
    
    # Get useful spots - ROI
    useful_spots <- get_ct_spots(slide, ct = ROI_ct, filter_prop = ROI_prop)
    
    if(length(useful_spots) < 25) {
      
      return(NULL)
      
    } 
    
    # Get predictor cell-types
    all_cts <- rownames(GetAssayData(slide, assay = "c2l"))
    all_cts <- all_cts[!(all_cts %in% c(ROI_ct, "prolif"))]
    
    # Get ct states that are useful
    
    # Here predictors and predicted values are considered separately
    # Easier to see the interaction between targets/predictors
    
    all_ct_states <- rownames(GetAssayData(slide, assay = "cell_states_pos"))
    
    useful_targets <- all_ct_states[all_ct_states %in% targets]
    
    # Define pipeline
    
    # Define assay of each view ---------------
    view_assays <- list("main" = target_assay,
                        "intra_pred" = "c2l",
                        "para_pred" = "c2l")
    
    # Define features of each view ------------
    view_features <- list("main" = useful_targets,
                          "intra_pred" = all_cts,
                          "para_pred" = all_cts)
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra",
                       "intra_pred" = "intra",
                       "para_pred" = "para")
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL,
                        "intra_pred" = NULL,
                        "para_pred" = 5)
    
    misty_out_alias <- paste0("./results/state_structure/", 
                              folder_label, 
                              "/mstate_")
    
    misty_out <- paste0(misty_out_alias, 
                        slide_id)
    
    run_misty_seurat(visium.slide = slide,
                     view.assays = view_assays,
                     view.features = view_features,
                     view.types = view_types,
                     view.params = view_params,
                     spot.ids = useful_spots,
                     out.alias = misty_out)
    
    misty_res_slide <- collect_results(misty_out)
    
    view_names <- misty_res_slide$importances$view %>% unique()
    
    plot_folder <- paste0(misty_out, "/plots")
    
    system(paste0("mkdir ", plot_folder))
    
    # Plot correlations for the sake of completion
    state_data_A <- GetAssayData(slide, assay = target_assay)[useful_targets, useful_spots, drop = F]
    state_data_B <- GetAssayData(slide, assay = "c2l")[all_cts, useful_spots]
    state_data <- rbind(state_data_A, state_data_B)
    
    cor_data_states <- cor(t(state_data)) %>%
      as.data.frame() %>%
      rownames_to_column("feature_a") %>%
      pivot_longer(-feature_a, names_to = "feature_b")
    
    write_csv(cor_data_states, file = paste0(plot_folder,"/state_cor.csv"))
    
    pdf(file = paste0(plot_folder, "/", slide_id, "_", "summary_plots.pdf"))
    
    mistyR::plot_improvement_stats(misty_res_slide)
    mistyR::plot_view_contributions(misty_res_slide)
    
    walk(view_names, function(v) {
      mistyR::plot_interaction_heatmap(misty_res_slide, v, cutoff = 0)
    })
    
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
  misty_outs <- misty_outs[grepl("mstate", misty_outs)]
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]])) %>%
    mutate(cor_res = paste0(misty_out_file, "/plots/state_cor.csv")) %>%
    mutate(cor_res = map(cor_res, read_csv))
  
  misty_res <- collect_results(misty_outs)
  
  model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                  grepl("multi.R2", measure)) %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]))
  
  importances_filtered <- misty_res$importances %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]])) %>%
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
  
  mean_importances %>% unnest() %>% write_csv(., paste0(misty_out_folder,"/summary_mean_importances.csv"))
  mean_cors %>% write_csv(., paste0(misty_out_folder,"/summary_mean_cors.csv"))
  
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


# Creates summary panels


plot_int_panels <- function(visium_slide, ROI_ct, param_list) {
  
  roi_spots <- get_ct_spots(visium_slide, 
                            ct = ROI_ct,
                            filter_prop = 0.1)
  
  visium_slide <- visium_slide[, roi_spots]
  
  param_list <- enframe(param_list,
                        name = "assay",
                        value = "features")
  
  plts <- pmap(param_list, function(assay, features) {
    
    DefaultAssay(visium_slide) <- assay
    
    all_plts <- SpatialFeaturePlot(visium_slide, features = features, combine = F,max.cutoff = "q99", min.cutoff = "q01")
    
    all_plts <- map(all_plts, function(plt) {
      
      plt + scale_fill_viridis()
      
    })
    
    all_plts <- plot_grid(plotlist = all_plts, align = "hv", ncol = 3)
    
  })
  
  final_panel <- plot_grid(plotlist = plts, ncol = 2, labels = c("predicted_states", "predictor_states", 
                                                                 "predictor_cts", "predictor_cts_para"),
                           align = "hv")
}

# Group correlations estimated from the pipeline

collect_correlations <- function(misty_out_folder, r2_filter, sample_dict) {
  
  misty_outs <- list.files(misty_out_folder, full.names = T)
  misty_outs <- misty_outs[grepl("mstate", misty_outs)]
  
  paramdf <- tibble(misty_out_file = misty_outs) %>%
    mutate(sample = strsplit(misty_out_file, 'mstate_') %>%
             map_chr(., ~.x[[2]])) %>%
    mutate(cor_res = paste0(misty_out_file, "/plots/state_cor.csv")) %>%
    mutate(cor_res = map(cor_res, read_csv))
  
  misty_res <- collect_results(misty_outs)
  
  model_performance <- misty_res$improvements %>% dplyr::filter(grepl("intra.R2", measure) | 
                                                                  grepl("multi.R2", measure)) %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]]))
  
  importances_filtered <- misty_res$importances %>%
    mutate(sample = strsplit(sample, 'mstate_') %>%
             map_chr(., ~.x[[2]])) %>%
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
  
  selected_cors <- best_performers %>%
    left_join(cors_filtered, 
              by = c("sample", 
                     "target" = "feature_a")) %>%
    na.omit() 
  
  return(selected_cors)
  
}

# Plot comparisons of correlations

plot_box_pw_cor <- function(cors, cor_comps, area = F) {
  
  cors_info <- cors %>%
    group_by(target, feature_b) %>%
    nest() %>%
    dplyr::rename("cor_data" = data)
  
  pw_info <- cor_comps %>%
    group_by(target, feature_b) %>%
    dplyr::select(group1, group2, p , p.adj) %>%
    nest() %>%
    dplyr::rename("pw_data" = data)
  
  boxplot_df <- left_join(cors_info, pw_info) %>%
    dplyr::filter(map_lgl(pw_data, ~ is.null(.x)) != TRUE) %>%
    mutate(bxplt = map2(cor_data, pw_data, function(crs, pw, area_flag = area) {
      
      # Here you filter areas that were tested
      if(area_flag) {
        groups_list <- unique(c(pw$group1, pw$group2))
        crs <- crs %>%
          dplyr::filter(major_labl %in% groups_list)
      }
      
      max_val <- max(crs$correlation) + 0.05
      
      if(area_flag) {
        
        box_plt <- ggplot(crs, aes(x = major_labl, y = correlation, color = major_labl))
        
      } else {
        
        box_plt <- ggplot(crs, aes(x = patient_group, y = correlation, color = patient_group))
        
      }
      
      box_plt <- box_plt +
        geom_boxplot() +
        geom_point() +
        ggpubr::stat_pvalue_manual(pw, label = "p.adj", 
                                   y.position = max_val, 
                                   step.increase = 0.1,
                                   tip.length = 0.01,size = 3) +
        theme_classic() +
        theme(axis.text = element_text(size = 11),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.5, 
                                         size = 11),
              legend.position = "none") +
        ylab("correlation") 
      
      
    }))
  
  return(boxplot_df)
}




