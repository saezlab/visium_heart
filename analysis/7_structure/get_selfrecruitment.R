# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' How do cells relate to each other?
#' Are there cells that tend to form communities?

library(tidyverse)
library(Seurat)
library(MISTy)
source("./visium/visiumtools/misty_utils.R")

# Get all data files ------------------------------------------------------------
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)
slide_files_full <- paste0(slide_files_folder, 
                           slide_files)
slide_ids <- gsub("[.]rds", "", slide_files)

# Putting cell2location scores and para/juxta representations in a single tibble ---
deconv_res <- tibble("slide_path" = slide_files_full,
                     "slide" = slide_ids) %>%
  dplyr::mutate(deconv_mats = map(slide_path, function(x) {
    
    slide_obj <- readRDS(x)
    DefaultAssay(slide_obj) = "c2l"
    
    tibble("ct_densities" = as.matrix(slide_obj[["c2l"]]@data) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type,
                          values_to = "absolute_score") %>%
             mutate(cell_type = gsub("-", "_", cell_type)),
           "ct_para" = get_para_matrix(visium_slide = slide_obj,
                                       para_assay = "c2l",
                                       para_features = rownames(slide_obj),
                                       l = 10) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type,
                          values_to = "para_score"),
           "ct_juxta" = get_juxta_matrix(visium_slide = slide_obj,
                                         juxta_assay = "c2l",
                                         juxta_features = rownames(slide_obj),
                                         neighbor_thr = 10) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type,
                          values_to = "juxta_score"))
    
    
  }))

deconv_res <- deconv_res %>%
  group_by(slide) %>%
  mutate(ct_data = map(deconv_mats, function(x){
    reduce(x, left_join, by = c("name","cell_type"))
  }))
  
deconv_res <- deconv_res %>%
  dplyr::select(slide,
                ct_data) %>%
  unnest()
  
# First generate an analysis ignoring the patients -----------

all_cors <- deconv_res %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(para_cor = map(data, function(x) {
    broom::tidy(cor.test(x$absolute_score, x$para_score))
  }),
  juxta_cor = map(data, function(x) {
    broom::tidy(cor.test(x$absolute_score, x$juxta_score))
  })
  ) %>%
  dplyr::select(-data) %>%
  dplyr::mutate(para_cor = map_dbl(para_cor, ~.x[["estimate"]]),
                juxta_cor = map_dbl(juxta_cor, ~.x[["estimate"]]))


pdf("./visium_results_manuscript/structure/self_relation_all.pdf",
    height = 5, width = 5)

ggplot(all_cors, 
       aes(x = para_cor,
           y = juxta_cor,
           label = cell_type)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  theme(axis.text = element_text(size = 11)) +
  ggtitle("Self relationship")

dev.off()

# Repeat it for each patient --------------------------------

all_cors_pat <- deconv_res %>%
  ungroup() %>%
  group_by(cell_type, slide) %>%
  nest() %>%
  mutate(para_cor = map(data, function(x) {
    broom::tidy(cor.test(x$absolute_score, x$para_score))
  }),
  juxta_cor = map(data, function(x) {
    broom::tidy(cor.test(x$absolute_score, x$juxta_score))
  })
  ) %>%
  dplyr::select(-data) %>%
  dplyr::mutate(para_cor = map_dbl(para_cor, ~.x[["estimate"]]),
                juxta_cor = map_dbl(juxta_cor, ~.x[["estimate"]])) 

# Annotate slides -------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(Visium, New.Ids) %>%
  dplyr::rename(slide = Visium,
                patient = New.Ids)

all_cors_pat <- all_cors_pat %>%
  left_join(condition_dictionary, by = "slide")

pdf("./visium_results_manuscript/structure/self_relation_patient.pdf",
    height = 8, width = 8)

print(ggplot(all_cors_pat, 
       aes(x = para_cor,
           y = juxta_cor,
           color = cell_type)) +
  geom_point() +
  theme(legend.position = "bottom") +
  facet_wrap(.~patient, nrow = 3))

           
dev.off()














