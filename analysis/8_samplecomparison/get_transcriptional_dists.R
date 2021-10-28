# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# This script generates distances between samples using pseudobulk profiles 
# of all cell-types

library(scater)
library(tidyverse)
library(rdist)
source("./analysis/utils/pseudobulk_utils.R")

all_mats <- readRDS("./processed_snrnaseq/integration/psxpat_integrated_rnasamples_filt.rds")
pat_anns <- readRDS("./markers/snrna_patient_anns_revisions.rds")[, c("patient_id", "patient_group")] %>%
  unique()

get_corrected_dist <- function(count_matrix, pat_anns) {
  
  norm_mat <- cpm_norm(count_matrix)
  
  meta_data <- tibble(patient_id = colnames(count_matrix)) %>%
    left_join(pat_anns) %>%
    arrange(patient_group) %>%
    group_by(patient_group) %>%
    nest() %>%
    dplyr::rename("anns" = data) %>%
    dplyr::mutate(nsamples = map_dbl(anns, nrow)) %>%
    dplyr::filter(nsamples > 1) %>%
    dplyr::mutate(gex = map(anns, function(ann) { # Subset expression
      cell_ids <- ann$patient_id
      mat <- norm_mat[, cell_ids]
    })) %>% # Calculate angular distance
    dplyr::mutate(dist_relationships = map(gex, get_ang_dist)) %>%
    dplyr::mutate(intra_mean_dist = map_dbl(dist_relationships, function(x) {
      mean(x$angulardist)
    }))
  
  # 1 group vs the other
  
  group_comps <- combn(meta_data$patient_group,2) %>%
    t() %>%
    as_tibble() %>%
    dplyr::rename("groupA" = V1, 
                  "groupB" = V2) %>%
    dplyr::mutate(std_factor = map2_dbl(groupA, groupB, function(a,b) {
      all_dists <- dplyr::filter(meta_data, patient_group %in% c(a,b)) %>%
        dplyr::select(dist_relationships) %>%
        unnest() %>%
        ungroup()
      
       mean(all_dists$angulardist)
       
    })) %>%
    dplyr::mutate(ang_distances = map2(groupA, groupB, function(a,b) {
      
      mats <- dplyr::filter(meta_data, patient_group %in% c(a,b)) %>%
        pull(gex)
      
      get_ang_dist_paired(gex_a = mats[[1]], gex_b = mats[[2]])
      
      })) %>% 
    unnest() %>%
    dplyr::mutate(std_angulardist = angulardist/std_factor)
  
  # 1 group vs the rest
  
  std_factor <- meta_data %>%
    ungroup() %>%
    dplyr::select(dist_relationships) %>%
    unnest() %>%
    summarise(mean(angulardist)) %>%
    pull()
  
  full_mat_distance <- get_ang_dist((reduce(meta_data$gex, cbind))) %>%
    dplyr::left_join(pat_anns, by = c("pat_a" = "patient_id")) %>%
    dplyr::left_join(pat_anns, by = c("pat_b" = "patient_id")) %>%
    dplyr::filter(patient_group.x != patient_group.y) %>%
    dplyr::mutate(std_angulardist = angulardist/std_factor)
  
  
  return(list("contrast" = group_comps, "general" = full_mat_distance))
  
}

get_ang_dist <- function(gex) {

  ang_dist <- pdist(t(gex), metric = "angular")
  ang_dist[upper.tri(ang_dist,diag = T)] <- NA
  ang_dist <- ang_dist %>%
    as.data.frame() %>%
    rownames_to_column("pat_a") %>%
    pivot_longer(-pat_a,
                 names_to = "pat_b", 
                 values_to = "angulardist") %>%
    na.omit()
  
} 


get_ang_dist_paired <- function(gex_a, gex_b) {
  
  ang_dist <- cdist(t(gex_a), t(gex_b), metric = "angular") %>%
    as.data.frame() %>%
    rownames_to_column("pat_a") %>%
    pivot_longer(-pat_a,
                 names_to = "pat_b", 
                 values_to = "angulardist") %>%
    na.omit()
  
  return(ang_dist)
  
}


# Main
trans_shifts <- all_mats %>%
  group_by(cell_type) %>%
  dplyr::mutate(angular_dists = map(count_matrix, get_corrected_dist, pat_anns = pat_anns)) %>%
  dplyr::select(angular_dists) %>%
  mutate(angular_dists = map(angular_dists, enframe)) %>%
  unnest()

# Do plots

# General
general_dat <- trans_shifts %>%
  dplyr::filter(name == "general") %>%
  unnest() 

cell_order <- general_dat %>%
  group_by(cell_type) %>%
  summarize(med_dist = median(std_angulardist)) %>%
  arrange(-med_dist) %>%
  pull(cell_type)

pdf("./results/sample_comparison/gen_transcriptional_shifts.pdf", height = 3, width = 4)
plot(general_dat %>%
  ggplot(aes(x = factor(cell_type,
                        levels = cell_order), 
             y = std_angulardist)) +
  geom_boxplot() +
  geom_hline(yintercept = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 11)) +
  xlab("") +
  ylab("transcriptional shift"))
dev.off()

# Contrasts
general_dat <- trans_shifts %>%
  dplyr::filter(name == "contrast") %>%
  unnest() %>%
  mutate(cell_type = factor(cell_type,
                            levels = cell_order),
         contrast_id = paste0(groupA,"_vs_",groupB))


pdf("./results/sample_comparison/contrast_transcriptional_shifts.pdf", height = 5, width = 4)
plot(general_dat %>%
       ggplot(aes(x = cell_type, 
                  y = std_angulardist)) +
       geom_boxplot() +
       geom_hline(yintercept = 1) +
       theme_classic() +
       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
             axis.text = element_text(size = 11)) +
       xlab("") +
       ylab("transcriptional shift") +
       facet_wrap(.~contrast_id,nrow = 3)
       ) 
dev.off()
  
