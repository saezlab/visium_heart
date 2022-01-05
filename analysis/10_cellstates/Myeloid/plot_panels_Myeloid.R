# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates panels of state to ct analysis
library(Seurat)
library(tidyverse)
library(cowplot)
library(viridis)
source("./analysis/utils/misty_pipeline.R")

# Define cell-type and get best in class summary

misty_out <- "./results/state_structure/Myeloid_ct/"

best_performers <- paste0(misty_out, "best_performers.csv") %>%
  read.csv()

bp_plt <- best_performers %>%
  ggplot(aes(x = target, y = sample, fill = value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf(paste0(misty_out, "best_performers.pdf"), width = 4, height = 5)

plot(bp_plt)

dev.off()

misty.results <- mistyR::collect_results(folders = paste0(misty_out,"mstate_", best_performers$sample %>% unique()))

pdf(paste0(misty_out, "view_contributions.pdf"), width = 4, height = 5)

mistyR::plot_view_contributions(misty.results)

dev.off()


# Get parameters

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/misty_red_objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("_mistyassays.rds", "", slide_files)

param_df <- tibble(visium_file = paste0(slide_files_folder,slide_files),
                   pdf_file = paste0(misty_out, "panels/", slide_ids, ".pdf"))

ROI_ct <- "Myeloid"

pwalk(param_df, function(visium_file, pdf_file) {
  print(visium_file)
  
  visium_slide <- readRDS(visium_file)
  
  cell_states <- GetAssayData(visium_slide, assay = "cell_states") %>%
    rownames()
  
  cell_states <- cell_states[grepl(ROI_ct, cell_states)]
  
  cell_states_predictor <- GetAssayData(visium_slide, assay = "cell_states_predictor") %>%
    rownames()
  
  cell_states_predictor <- cell_states_predictor[grepl(ROI_ct, cell_states_predictor)]
  
  cell_states <- cell_states[grepl(ROI_ct, cell_states)]
  
  c2l <- GetAssayData(visium_slide, assay = "c2l") %>%
    rownames()
  
  c2l <- c2l[!c2l %in% c(ROI_ct, "prolif")]
  
  c2l_para_predictor <- GetAssayData(visium_slide, assay = "c2l_para_predictor") %>%
    rownames()
  
  c2l_para_predictor <- c2l_para_predictor[!c2l_para_predictor %in% c(ROI_ct, "prolif")]
  
  param_list <- list("cell_states" = cell_states,
                     "cell_states_predictor" = cell_states_predictor,
                     "c2l" = c2l,
                     "c2l_para_predictor" = c2l_para_predictor)
  
  
  panels <- plot_int_panels(visium_slide = visium_slide,ROI_ct = ROI_ct, param_list = param_list)
  
  pdf(pdf_file, width = 20, height = 16)
  
  plot(panels)
  
  dev.off()
  
  
})

# Repeat for Li's annotation


misty_out <- "./results/state_structure/Myeloid_ct_ms/"

best_performers <- paste0(misty_out, "best_performers.csv") %>%
  read.csv()

bp_plt <- best_performers %>%
  ggplot(aes(x = target, y = sample, fill = value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf(paste0(misty_out, "best_performers.pdf"), width = 4, height = 5)

plot(bp_plt)

dev.off()

# Get parameters

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds")
slide_files_folder <- "./processed_visium/misty_red_objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("_mistyassays.rds", "", slide_files)

param_df <- tibble(visium_file = paste0(slide_files_folder,slide_files),
                   pdf_file = paste0(misty_out, "panels/", slide_ids, ".pdf"))

ROI_ct <- "Myeloid"

pwalk(param_df, function(visium_file, pdf_file) {
  print(visium_file)
  
  visium_slide <- readRDS(visium_file)
  
  cell_states <- GetAssayData(visium_slide, assay = "cell_states") %>%
    rownames()
  
  cell_states <- cell_states[grepl(ROI_ct, cell_states)]
  
  cell_states_predictor <- GetAssayData(visium_slide, assay = "cell_states_predictor") %>%
    rownames()
  
  cell_states_predictor <- cell_states_predictor[grepl(ROI_ct, cell_states_predictor)]
  
  cell_states <- cell_states[grepl(ROI_ct, cell_states)]
  
  c2l <- GetAssayData(visium_slide, assay = "c2l") %>%
    rownames()
  
  c2l <- c2l[!c2l %in% c(ROI_ct, "prolif")]
  
  c2l_para_predictor <- GetAssayData(visium_slide, assay = "c2l_para_predictor") %>%
    rownames()
  
  c2l_para_predictor <- c2l_para_predictor[!c2l_para_predictor %in% c(ROI_ct, "prolif")]
  
  param_list <- list("cell_states_ms_target" = cell_states,
                     "cell_states_ms_predictor" = cell_states_predictor,
                     "c2l" = c2l,
                     "c2l_para_predictor" = c2l_para_predictor)
  
  
  panels <- plot_int_panels(visium_slide = visium_slide,ROI_ct = ROI_ct, param_list = param_list)
  
  pdf(pdf_file, width = 20, height = 16)
  
  plot(panels)
  
  dev.off()
  
  
})











