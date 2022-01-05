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
predictors <- c("Fib-Fib-0", "Fib-Fib-3",
                 "Fib-Fib-SCARA5","Fib-Myofib")

target_list <- c("Myeloid-DCs-FLT3-ITGAX", "Myeloid-LYVE-FOLR-Macrophages",
                "Myeloid-LYVE-PLTP-Macrophages", "Myeloid-Monocyte-CCL18",
                "Myeloid-Monocyte-SPP1")

# Myeloid ROI -----------------------------------------------

# Joint analysis ------------------------

run_state_ppline(ROI_ct = "Myeloid",
                 ROI_prop = 0.1,
                 mask_by_prop = F,
                 mask_threshold = 0.1,
                 folder_label = "Myeloid_Fib",
                 targets =  target_list,
                 predictors = predictors,
                 join_trgts_pred = F,
                 target_assay = "cell_states",
                 state_origin = "cell_states")#Defines region of interest,

misty_out_folder <- "./results/state_structure/Myeloid_Fib/"
performance_all_misty(misty_out_folder, r2_filter = 0.1)

# Weighted mapping by proportions

run_state_ppline(ROI_ct = "Myeloid",
                 ROI_prop = 0.1,
                 mask_by_prop = F, # Weighted by cell proportions
                 mask_threshold = 0.1,
                 folder_label = "Myeloid_Fib_pos",
                 targets =  target_list,
                 predictors = predictors,
                 join_trgts_pred = F,
                 target_assay = "cell_states",
                 state_origin = "cell_states")#Defines region of interest





