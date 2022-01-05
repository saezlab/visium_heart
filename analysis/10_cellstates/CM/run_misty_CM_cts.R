# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Test a simplified version of the spatial analysis of interacting cells of interest

library(tidyverse)
library(Seurat)
library(mistyR)
source("./analysis/utils/misty_utilities.R")
source("./analysis/utils/misty_pipeline.R")

future::plan(future::multisession)

# Main ------------------------------------------------------------------------
# Getting sample annotations --------------------------------------------------
sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

# Targets
target_list <- c("CM-damaged-CM", "CM-healthy-CM", "CM-vCM-4")

# CM ROI -----------------------------------------------

# This is wmean in ROI of interest
run_state_ppline_ct(ROI_ct = "CM",
                    ROI_prop = 0.1,
                    mask_by_prop = F,
                    mask_threshold = 0.1,
                    folder_label = "CM_ct",
                    targets =  target_list,
                    target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/CM_ct/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)

# This is module_score in ROI of interest
run_state_ppline_ct(ROI_ct = "CM",
                    ROI_prop = 0.1,
                    mask_by_prop = T,
                    mask_threshold = NULL,
                    folder_label = "CM_ct_ms",
                    targets =  target_list,
                    target_assay = "cell_states_pos",
                    state_origin = "cell_states_ms")

misty_out_folder <- "./results/state_structure/CM_ct_ms/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)

# Simpler non-redundant object

# Targets
target_list <- c("CM-damaged-CM") 

# CM ROI -----------------------------------------------

# This is wmean in ROI of interest
run_state_ppline_ct_simple(ROI_ct = "CM",
                    ROI_prop = 0.1,
                    mask_by_prop = F,
                    mask_threshold = 0.1,
                    folder_label = "CM_ct_damaged",
                    targets =  target_list,
                    target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/CM_ct_damaged/"
performance_all_misty(misty_out_folder, r2_filter = 20)

# CM ROI -----------------------------------------------
target_list <- c("CM-healthy-CM")

# This is wmean in ROI of interest
run_state_ppline_ct_simple(ROI_ct = "CM",
                           ROI_prop = 0.1,
                           mask_by_prop = F,
                           mask_threshold = 0.1,
                           folder_label = "CM_ct_healthy",
                           targets =  target_list,
                           target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/CM_ct_healthy/"
performance_all_misty(misty_out_folder, r2_filter = 20)





# This is module_score in ROI of interest
run_state_ppline_ct_simple(ROI_ct = "CM",
                    ROI_prop = 0.1,
                    mask_by_prop = T,
                    mask_threshold = NULL,
                    folder_label = "CM_ct_ms_simple",
                    targets =  target_list,
                    target_assay = "cell_states_pos",
                    state_origin = "cell_states_ms")

misty_out_folder <- "./results/state_structure/CM_ct_ms_simple/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)


# run_state_ppline_ct(ROI_ct = "CM",
#                     ROI_prop = 0.1,
#                     mask_by_prop = F,
#                     mask_threshold = 0.1,
#                     folder_label = "CM_ct_pos",
#                     targets =  target_list,
#                     target_assay = "cell_states_pos")
# 
# misty_out_folder <- "./results/state_structure/CM_ct_pos/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)

# Examining importances

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



