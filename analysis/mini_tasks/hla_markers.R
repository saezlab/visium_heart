# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create hmaps per patient, per cell-type of interest
library(Seurat)
library(scater)
library(ComplexHeatmap)
library(tidyverse)
source("./analysis/utils/pseudobulk_utils.R")

scell_data <- readRDS("./visium_results_manuscript/ct_data/fibroblasts_states.rds")

scell_data_ps <- sumCountsAcrossCells(x = as.matrix(scell_data@assays[["RNA"]]@counts),
                     ids = scell_data@meta.data[, "orig.ident"])

scell_data_ps <- assay(scell_data_ps)

scell_data_ps <- scell_data_ps[grepl("^HLA", rownames(scell_data_ps)), ]

condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(snRNA, New.Ids) %>%
  dplyr::rename(orig.ident = snRNA,
                patient = New.Ids) %>%
  write.table("./visium_results_manuscript/mini_tasks/condition_dictionary.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)

write.table(scell_data_ps, file = "./visium_results_manuscript/mini_tasks/HLA_pseudobulk_patients.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

hmap_plt <- ComplexHeatmap::Heatmap(scell_data_ps,
                                    name = "expr",
                                    show_row_dend = F)
