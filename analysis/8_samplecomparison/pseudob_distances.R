# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# This script generates distances between samples using pseudobulk profiles 
# of all cell-types

library(scater)
library(tidyverse)
library(philentropy)
library(cowplot)
source("./analysis/utils/pseudobulk_utils.R")

pseudobulk_data <- readRDS("./visium_results_manuscript/integration/ps_patients_ct.rds")[[1]][["gex"]]
atlas_meta <- readRDS("./visium_results_manuscript/integration/ps_patients_ct.rds")[[1]][["annotations"]]

# Get sample information ----------------------------------------------------------------------

# ID names
sample_dict <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(snRNA, New.Ids) %>%
  dplyr::rename(orig.ident = snRNA,
                patient_sample = New.Ids) %>%
  dplyr::mutate(patient = map_chr(strsplit(patient_sample, "_"), 
                                  ~.x[[1]]),
                area = map_chr(strsplit(patient_sample, "_"), 
                               ~.x[[2]]))

# Cell-type composition
cell_props <- atlas_meta %>%
  group_by(orig.ident, cell_type) %>%
  summarize(ncells = length(cell_type)) %>%
  dplyr::mutate(all_sample_cells = sum(ncells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cell_prop = ncells/all_sample_cells) %>%
  dplyr::select(orig.ident, cell_type, cell_prop)


# Get matrices per cell_type ----------------------------------------------------------------------
col_meta <- colData(pseudobulk_data) %>%
  as.data.frame()

cts <- set_names(unique(col_meta %>% pull(cell_type)))

cell_type_divergences <- map(cts, function(ct) {
  
  # First, get indexes
  ct_ix <- which(col_meta %>% pull(cell_type) == ct)
  
  ct_dat <- assay(pseudobulk_data)[, ct_ix]
  colnames(ct_dat) <- (col_meta %>% pull(orig.ident))[ct_ix]
  
  ct_dat <- edgeR_filtering(ct_dat, min.count = 0.2)
  
}) %>% 
  enframe(name = "cell_type", value = "count_matrix")

# Calculate JSD distances from counts -----------------------------------

JSD_mod <- function(cmat) {
  
  ids <- colnames(cmat)
  
  jsd_dist <- as.data.frame(philentropy::JSD(t(cmat), 
                                             est.prob = "empirical"))
  
  rownames(jsd_dist) <- colnames(jsd_dist) <- ids
  
  # Make it distance ----------------------------------
  jsd_dist <- jsd_dist ** (1/2)

  jsd_dist <- as.data.frame(jsd_dist) %>%
    rownames_to_column("SampleA") %>%
    pivot_longer(-SampleA, 
                 names_to = "SampleB",
                 values_to = "jsd") %>%
    na.omit()
  
  return(jsd_dist)
  
}

# Calculate divergence of each cell-type for each pair of samples
cell_type_divergences <- cell_type_divergences %>% 
  group_by(cell_type) %>%
  dplyr::mutate(jsd_dist = map(count_matrix, JSD_mod)) %>%
  dplyr::select(cell_type, jsd_dist) %>%
  unnest() %>%
  ungroup()

# Weight divergences by proportions
cell_type_divergences <- left_join(cell_type_divergences, cell_props,
                                   by = c("SampleA" = "orig.ident",
                                          "cell_type")) %>%
  dplyr::rename(cell_prop_A = cell_prop) %>%
  left_join(cell_props,
            by = c("SampleB" = "orig.ident",
                   "cell_type")) %>%
  dplyr::rename(cell_prop_B = cell_prop) %>%
  dplyr::mutate(min_prop = ifelse(cell_prop_A <= cell_prop_B, 
                                  cell_prop_A, cell_prop_B))

# Sample divergences ---------------------------------------
sample_divergences <- cell_type_divergences %>%
  dplyr::mutate(weighted_jsd = min_prop * jsd) %>%
  group_by(SampleA, SampleB) %>%
  summarize(patient_jsd = sum(weighted_jsd)) %>%
  pivot_wider(names_from = SampleA,
              values_from = patient_jsd)

# Multidimensional scaling using JSD distances
get_JSD_mds <- function(divergences) {
  rows <- divergences[[1]]
  sample_divergences <- as.matrix(divergences[, -1])
  rownames(sample_divergences) <- rows
  
  mds_fit <- as.data.frame(cmdscale(as.dist(sample_divergences), 
                                    eig=TRUE, k=2)$points) %>%
    rownames_to_column("orig.ident") %>% 
    left_join(sample_dict)
  
  return(mds_fit)
}

# MDS of complete profile -----------------------------------------

mds_samples_all <- get_JSD_mds(sample_divergences)

all_jsd_plt <- ggplot(mds_samples_all, aes(x = V1,
                                           y = V2,
                                           color = area,
                                           label = patient)) + 
  geom_point(size = 2) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  ggtitle("Weighted JSD (all profile)")

# Do it per cell-type

cell_type_divergences <- cell_type_divergences %>%
  dplyr::select(cell_type, SampleA, SampleB, jsd) %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(cell_type_mds = map(data, function(x) {
    get_JSD_mds(pivot_wider(x, names_from = SampleA,
                            values_from = jsd))
  })) %>%
  dplyr::mutate(mds_plts = map2(cell_type_mds, cell_type, function(x, y) {
    ggplot(x, aes(x = V1, 
                  y = V2, 
                  color = area,
                  label = patient)) + 
      geom_point(size = 1) +
      ggrepel::geom_text_repel() +
      theme_classic() +
      ggtitle(paste0("JSD ", y))
  }))

# generate panel

JSD_plts = cowplot::plot_grid(all_jsd_plt, cell_type_divergences$mds_plts[[1]], cell_type_divergences$mds_plts[[2]],
                     cell_type_divergences$mds_plts[[3]], cell_type_divergences$mds_plts[[4]], cell_type_divergences$mds_plts[[5]],
                     cell_type_divergences$mds_plts[[6]], cell_type_divergences$mds_plts[[7]], cell_type_divergences$mds_plts[[8]],
                     cell_type_divergences$mds_plts[[9]],
                     nrow = 3, ncol = 4, align = "hv")

pdf(file = "visium_results_manuscript/sample_comparison/JSD_distances.pdf", width = 17, height = 11)

plot(JSD_plts)

dev.off()







