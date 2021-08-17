# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we evaluate the deconvolution runs 
#' In theory one already has the counts in the slide objects after
#' running single_nmf

library(tidyverse)
library(Seurat)

# Get all data files
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)
slide_files_full <- paste0(slide_files_folder, 
                           slide_files)
slide_ids <- gsub("[.]rds", "", slide_files)

# Patient annotation -------------------------------------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(Visium, New.Ids) %>%
  dplyr::rename(slide = Visium,
                patient = New.Ids)

# Tibble of all data
# c2l: cell2location
# spotlight: spotlight

deconv_res <- tibble("slide_path" = slide_files_full,
                     "slide" = slide_ids) %>%
  dplyr::mutate(deconv_mats = map(slide_path, function(x) {
    
    slide_obj <- readRDS(x)
    
    tibble("cell2location" = as.matrix(slide_obj[["c2l"]]@data) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type),
           "spotlight" = as.matrix(slide_obj[["spotlight"]]@data) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type))
    
    
  }))

# Here I keep in long format c2l and spotlight runs for further filtering

deconv_res <- deconv_res %>% 
  dplyr::select(slide, deconv_mats) %>%
  dplyr::mutate(spotlight = map(deconv_mats, ~ .x[["spotlight"]]),
                cell2location = map(deconv_mats, ~ .x[["cell2location"]])) %>%
  dplyr::select(slide, spotlight, cell2location) %>%
  dplyr::left_join(condition_dictionary)

# Spotlight inspection

spoligth_deconv <- deconv_res %>%
  dplyr::select(patient, spotlight) %>%
  unnest() %>%
  dplyr::filter(value >= 0.1) %>%
  ggplot(aes(y = cell_type,
             x = value)) +
  geom_violin() +
  xlim(c(0,1)) +
  facet_wrap(.~patient)

print(spoligth_deconv)

# Cell2location

c2l_deconv <- deconv_res %>%
  dplyr::select(patient, cell2location) %>%
  unnest() %>%
  ggplot(aes(y = cell_type,
             x = value)) +
  geom_violin() +
  facet_wrap(.~patient)

print(c2l_deconv)

# Transform them to proportions?
# Same filtering as spotlight
c2l_deconv_prop <- deconv_res %>%
  dplyr::select(patient, cell2location) %>%
  unnest() %>%
  group_by(patient, name) %>%
  dplyr::mutate(prop = value/(sum(value))) %>%
  ungroup() %>%
  dplyr::filter(prop >= 0.1) %>%
  ggplot(aes(y = cell_type,
             x = prop)) +
  geom_violin() +
  xlim(c(0,1)) +
  facet_wrap(.~patient)












