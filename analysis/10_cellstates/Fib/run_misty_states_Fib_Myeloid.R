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
target_list <- c("Fib-Fib-0", "Fib-Fib-3",
                 "Fib-Fib-SCARA5","Fib-Myofib")

predictors <- c("Myeloid-DCs-FLT3-ITGAX", "Myeloid-LYVE-FOLR-Macrophages",
                "Myeloid-LYVE-PLTP-Macrophages", "Myeloid-Monocyte-CCL18",
                "Myeloid-Monocyte-SPP1")

# Main
# Fibrosis ROI -----------------------------------------------
# Target scores aren't modified
run_state_ppline(ROI_ct = "Fib",
                 ROI_prop = 0.1,
                 mask_by_prop = F,
                 mask_threshold = 0.1,
                 folder_label = "Fib_Myeloid",
                 targets =  target_list,
                 predictors = predictors,
                 join_trgts_pred = F,
                 target_assay = "cell_states",
                 state_origin = "cell_states")#Defines region of interest,

misty_out_folder <- "./results/state_structure/Fib_Myeloid/"
performance_all_misty(misty_out_folder, r2_filter = 10)

# # Targets are positive and masked
# run_state_ppline(ROI_ct = "Fib",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Fib_Myeloid_pos",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = F,
#                  target_assay = "cell_states_pos",
#                  state_origin = "cell_states")#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Fib_Myeloid_pos/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)

# 
# # Weighted mapping by proportions
# 
# run_state_ppline(ROI_ct = "Fib",
#                  ROI_prop = 0.1,
#                  mask_by_prop = T, # Weighted by cell proportions
#                  mask_threshold = NULL,
#                  folder_label = "Fib_Myeloid_prop",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = T)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Fib_Myeloid_prop/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# # Separated analysis
# 
# run_state_ppline(ROI_ct = "Fib",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Fib_Myeloid_sep",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Fib_Myeloid_sep/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# run_state_ppline(ROI_ct = "Fib",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Fib_Myeloid_sep_unmasked",
#                  targets =  target_list,
#                  predictors = predictors,
#                  target_assay = "cell_states",
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Fib_Myeloid_sep_unmasked/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# run_state_ppline(ROI_ct = "Fib",
#                  ROI_prop = 0.1,
#                  mask_by_prop = T,
#                  mask_threshold = NULL,
#                  folder_label = "Fib_Myeloid_sep_prop",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Fib_Myeloid_sep_prop/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# # Myeloid ROI -----------------------------------------------
# 
# # Joint analysis ------------------------
# 
# run_state_ppline(ROI_ct = "Myeloid",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Myeloid_Fib",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = T)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# # Weighted mapping by proportions
# 
# run_state_ppline(ROI_ct = "Myeloid",
#                  ROI_prop = 0.1,
#                  mask_by_prop = T, # Weighted by cell proportions
#                  mask_threshold = NULL,
#                  folder_label = "Myeloid_Fib_prop",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = T)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib_prop/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# # Separated analysis
# 
# run_state_ppline(ROI_ct = "Myeloid",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Myeloid_Fib_sep",
#                  targets =  predictors,
#                  predictors = target_list,
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib_sep/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# run_state_ppline(ROI_ct = "Myeloid",
#                  ROI_prop = 0.1,
#                  mask_by_prop = F,
#                  mask_threshold = 0.1,
#                  folder_label = "Myeloid_Fib_sep_unmasked",
#                  targets =  predictors,
#                  predictors = target_list,
#                  target_assay = "cell_states",
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib_sep_unmasked/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# 
# 
# run_state_ppline(ROI_ct = "Myeloid",
#                  ROI_prop = 0.1,
#                  mask_by_prop = T,
#                  mask_threshold = NULL,
#                  folder_label = "Myeloid_Fib_sep_prop",
#                  targets =  target_list,
#                  predictors = predictors,
#                  join_trgts_pred = F)#Defines region of interest,
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib_sep_prop/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# 
# best_performers <- read_csv("./results/state_structure/Fib_Myeloid/best_performers.csv") %>%
#   left_join(sample_dict, by = c("sample" = "sample_id"))
# 
# best_performers %>%
#   ggplot(aes(x = patient_group, y = value)) +
#   geom_boxplot() +
#   facet_wrap(.~target)
# 
# best_performers %>%
#   group_by(target, patient_group) %>%
#   summarise(n())
# 
# # Be more stringent
# 
# run_state_ppline(ct = "Myeloid", #Defines region of interest
#                  folder_label = "Myeloid_Fib", 
#                  model_states = c(predictors,target_list))
# 
# misty_out_folder <- "./results/state_structure/Myeloid_Fib/"
# performance_all_misty(misty_out_folder, r2_filter = 0.1)
# 
# best_performers <- read_csv("./results/state_structure/Fib_Myeloid/best_performers.csv") %>%
#   left_join(sample_dict, by = c("sample" = "sample_id"))
# 
