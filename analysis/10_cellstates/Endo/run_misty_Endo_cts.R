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
target_list <- c("Endo-Arterial-Endo", "Endo-Capillary-Endo", "Endo-Endocardial-Endo",
                 "Endo-Lymphatic-Endo", "Endo-Venous-Endo")

# CM ROI -----------------------------------------------

run_state_ppline_ct(ROI_ct = "Endo",
                    ROI_prop = 0.1,
                    mask_by_prop = F,
                    mask_threshold = 0.1,
                    folder_label = "Endo_ct",
                    targets =  target_list,
                    target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/Endo_ct/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)


run_state_ppline_ct(ROI_ct = "Endo",
                    ROI_prop = 0.1,
                    mask_by_prop = F,
                    mask_threshold = 0.1,
                    folder_label = "Endo_ct_ms",
                    targets =  target_list,
                    target_assay = "cell_states_pos",
                    state_origin = "cell_states_ms")

misty_out_folder <- "./results/state_structure/Endo_ct_ms/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)


# run_state_ppline_ct(ROI_ct = "Endo",
#                     ROI_prop = 0.1,
#                     mask_by_prop = F,
#                     mask_threshold = 0.1,
#                     folder_label = "Endo_ct_pos",
#                     targets =  target_list,
#                     target_assay = "cell_states_pos")
# 
# misty_out_folder <- "./results/state_structure/Endo_ct_pos/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)

