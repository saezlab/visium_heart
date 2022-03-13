# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we estimate cell compositions from visium

library(tidyverse)
library(Seurat)

# Get all data files -----------------------------------------------------------------------
slide_files_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"
slide_files <- list.files(slide_files_folder)
slide_files_full <- paste0(slide_files_folder, 
                           slide_files)
slide_ids <- gsub("[.]rds", "", slide_files)

# Patient annotation -------------------------------------------------------------------------
visium_anns <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  dplyr::rename("visium_sample_id" = sample_id)

deconv_res <- tibble("slide_path" = slide_files_full,
                     "visium_sample_id" = slide_ids) %>%
  dplyr::mutate(deconv_mats = map(slide_path, ~ readRDS(.x) %>%
                                    as.data.frame() %>%
                                    rownames_to_column("spot_id") %>%
                                    pivot_longer(-spot_id, 
                                                 names_to = "cell_type",
                                                 values_to = "c2l_value"))) %>%
  dplyr::select(-slide_path) %>%
  unnest()

# Then we filter all location scores that represent less than 1 cell per spot given our priors
# prior from nuclei quantification

deconv_res <- deconv_res %>%
  group_by(visium_sample_id, spot_id) %>%
  mutate(n_cells_spot =  sum(c2l_value)) %>%
  mutate(c2l_value_prop = c2l_value / n_cells_spot) %>%
  dplyr::filter(n_cells_spot > 0) %>%
  ungroup()

# Here we add the patient information

deconv_res <- left_join(deconv_res, visium_anns)

spatial_props <- deconv_res %>%
  group_by(patient_id, cell_type) %>%
  summarize(sp_n_cells = sum(c2l_value)) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  mutate(sp_prop_cells = sp_n_cells/sum(sp_n_cells))

write.table(spatial_props, file = "./results/compositions/spatial_compositions.txt", col.names = T, row.names = F, quote = F, sep = "\t")











