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
target_list <- c("Myeloid-DCs-FLT3-ITGAX", "Myeloid-LYVE-FOLR-Macrophages",
                 "Myeloid-LYVE-PLTP-Macrophages", "Myeloid-Monocyte-CCL18",
                 "Myeloid-Monocyte-SPP1")

# Fibrosis ROI -----------------------------------------------

run_state_ppline_ct(ROI_ct = "Myeloid",
                    ROI_prop = 0.1,
                    mask_by_prop = F,
                    mask_threshold = 0.1,
                    folder_label = "Myeloid_ct",
                    targets =  target_list,
                    target_assay = "cell_states")

misty_out_folder <- "./results/state_structure/Myeloid_ct/"
performance_all_misty(misty_out_folder, r2_filter = 25)

# run_state_ppline_ct(ROI_ct = "Myeloid",
#                     ROI_prop = 0.1,
#                     mask_by_prop = F,
#                     mask_threshold = 0.1,
#                     folder_label = "Myeloid_ct_pos",
#                     targets =  target_list,
#                     target_assay = "cell_states_pos")
# 
# misty_out_folder <- "./results/state_structure/Myeloid_ct_pos/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)



