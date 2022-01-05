# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Get joint annotation of cell-states

library(Seurat)
library(tidyverse)

annotation_list <- readRDS("./processed_snrnaseq/cell_states/cellstate_annotation_list.rds")

all_annotations <- enframe(annotation_list) %>% 
  unnest() %>%
  dplyr::select(-name)

sc_data <- readRDS("./processed_snrnaseq/integration/integrated_rnasamples_ann.rds")

cell_state_annotation <- sc_data@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("raw_id") %>%
  left_join(all_annotations) %>%
  dplyr::filter(!is.na(annotation))

sc_data <- sc_data[, cell_state_annotation$raw_id]

sc_data$annotation <- cell_state_annotation$annotation

saveRDS(sc_data, file = "./processed_snrnaseq/cell_states/integrated_rnasamples_ann_wstates.rds")