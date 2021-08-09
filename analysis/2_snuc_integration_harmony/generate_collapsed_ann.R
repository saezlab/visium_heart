# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create annotation of states after filtering/collapsing

library(tidyverse)

# This is the original annotation after cell-state identification
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_wstates.rds")[[1]]
# Adding vsmcs
scell_meta <- mi_atlas[["annotations"]] %>%
  rownames_to_column("cell_id") %>% 
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id, split = "-"), ~.x[1])) %>%
  dplyr::mutate(cell_type = ifelse(deconv_col == "VSMCs",
                                   "VSMCs", cell_type)) %>%
  dplyr::select(cell_id, orig.ident, cell_type, deconv_col)

# This is the filtering performed in reduce_states.R
collapsing_stats <- readRDS("./visium_results_manuscript/ct_data/collapsing_data.rds")

patient_prop_filtering <- collapsing_stats$patient_filter

state_corr_filtering <- collapsing_stats$correlation_filter

# First tag all cell-states that were filtered because of patient specificity ------------------

scell_meta <- scell_meta %>%
  mutate(filtered_by_patspc = FALSE) %>%
  mutate(filtered_by_patspc = ifelse( (cell_type %in% unique(patient_prop_filtering$cell_type)) &
                                      !(deconv_col %in% unique(patient_prop_filtering$deconv_col)),
                                      TRUE, FALSE))

# Then filter cells that look like everything else (undefined) ---------------------------------

undefined_states <- state_corr_filtering %>% 
  dplyr::select(-simil_mat_plt) %>%
  unnest() %>%
  dplyr::mutate(n_states = length(unique(cell_state_a))) %>%
  group_by(cell_state_a) %>%
  dplyr::mutate(n_simil = sum(identical_states)) %>%
  ungroup() %>% 
  dplyr::select(cell_state_a, n_states, n_simil) %>%
  unique() %>%
  dplyr::filter(n_simil == (n_states - 1)) %>%
  pull(cell_state_a)

scell_meta <- scell_meta %>%
  dplyr::mutate(filtered_by_undefinition = ifelse(deconv_col %in% undefined_states, 
                                                  TRUE, FALSE))

# Create a new mapping for cells that need modification

new_mapping <- read.table("./markers/collapsed_states_ann.txt", 
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F)

scell_meta <- left_join(scell_meta, new_mapping,
                        by = c("deconv_col" = "state")) %>%
  dplyr::mutate(deconv_col = ifelse(is.na(grouped_state), 
                                    deconv_col, grouped_state)) %>%
  dplyr::select(-cell_type)


write.table(scell_meta, sep = "\t", quote = F,
            row.names = F, col.names = T,
            file = "./visium_results_manuscript/ct_data/states_annotations_aftercollapse.txt")


