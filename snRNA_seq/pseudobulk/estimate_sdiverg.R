# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Compares pseudobulk profiles of major cell_types

library(scater)
library(tidyverse)
library(rdist)
library(philentropy)
library(ggrepel)
library(cowplot)

# Read raw data -----------------------------------------------------------------------
pseudobulk_data <- readRDS("visium_results_manuscript/pseudobulk/mi_pseudobulk.rds")
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T)

# Get additional annotations ---------------------------------------------------------
condition_dictionary <- condition_dictionary %>%
  dplyr::select(snRNA, New.Ids) %>%
  dplyr::mutate(patient = map_chr(strsplit(New.Ids, "_"), ~ .x[1])) %>%
  dplyr::mutate(zone = map_chr(strsplit(New.Ids, "_"), ~ .x[2])) %>%
  dplyr::select(-New.Ids) %>%
  arrange(zone) %>%
  dplyr::mutate("patient_id" = paste0(patient, "-",zone))

# Fix names from pseudobulk ------------------------------------------------------
names(pseudobulk_data) <- map_chr(strsplit(map_chr(strsplit(names(pseudobulk_data), 
                                                           "/"),
                                                  dplyr::last), 
                                          "[.]"), 
                                 ~ .x[1])

# Get pseudobulk profiles of grouping variable of interest -----------------------
pseudobulk_data <- map(pseudobulk_data, ~ .x[["major_cell_type"]])

# Get counts from selected data -----------------------
cell_counts <- map(pseudobulk_data, function(x) {
  tibble("cell_id" = x$ids, 
         "counts" = x$ncells,
         "proportions" = x$ncells/sum(x$ncells))
}) %>% enframe(name = "snRNA") %>%
  unnest() %>%
  left_join(condition_dictionary, by = "snRNA") %>%
  dplyr::mutate("patient_id" = paste0(patient, "-",zone))

# Visualization of major cell types --------------------------
prop_overview <- ggplot(cell_counts, 
       aes(fill = cell_id, 
           x = counts, 
           y = factor(patient_id,
                      levels = condition_dictionary$patient_id))) + 
  geom_bar(position="fill", 
           stat="identity") +
  theme(legend.position = "bottom") +
  ylab("") + xlab("proportions")

# Get pseudobulk profiles of grouping variable of interest -----------------------
pseudobulk_data <- map(pseudobulk_data, assay)
scale_factor <- 10000
pseudobulk_data_normalized <- map(pseudobulk_data, function(x) {
  
  norm_mat <- log1p((x + 1)/sum(x) * scale_factor)
  
  return(norm_mat)
  
})

# Tidying: we need matrices per cell_type ----------------------

pseudobulk_data_long <- map(pseudobulk_data, function(x) {
  
    long_matrix <- x %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>%
    pivot_longer(-gene, 
                 names_to = "cell_id",
                 values_to = "counts")
    
}) %>% enframe(name = "snRNA") %>%
  unnest()

# Define function to make matrices all over again -------------------
fromtibble_togex <- function(data) {
  
  wide_form <- data %>%
    pivot_wider(names_from = snRNA,
                values_from = counts)
  
  expr_mat <- as.matrix(wide_form[,-1])
  
  rownames(expr_mat) <- wide_form$gene
  
  return(na.omit(expr_mat))
  
}
# Define function to wrap Jensen-Shannon divergence
# to deal with dependencies (distance relative to other points)
JSD_mod <- function(cmat) {
  
  ids <- colnames(cmat)
  
  jsd_dist <- as.data.frame(philentropy::JSD(t(cmat), 
                                             est.prob = "empirical"))
  
  rownames(jsd_dist) <- colnames(jsd_dist) <- ids
  
  #jsd_dist[!upper.tri(jsd_dist)] <- NA
  
  jsd_dist <- as.data.frame(jsd_dist) %>%
    rownames_to_column("SampleA") %>%
    pivot_longer(-SampleA, 
                 names_to = "SampleB",
                 values_to = "jsd") %>%
    na.omit()
  
  return(jsd_dist)
  
}

# Calculate divergence of each cell-type for each pair of samples 

cell_type_divergences <- pseudobulk_data_long %>%
  group_by(cell_id ) %>%
  nest() %>%
  dplyr::mutate(count_matrix = map(data, fromtibble_togex)) %>%
  dplyr::select(cell_id , count_matrix) %>%
  dplyr::mutate(npatients = map_int(count_matrix, ncol)) %>%
  dplyr::filter(npatients > 1) %>%
  dplyr::mutate(jsd_dist = map(count_matrix, JSD_mod)) %>%
  dplyr::select(cell_id, jsd_dist) %>%
  unnest() %>%
  ungroup()

# Weight divergences by proportions
cell_type_divergences <- left_join(cell_type_divergences, cell_counts,
                                   by = c("SampleA" = "snRNA",
                                          "cell_id")) %>%
  dplyr::rename(cell_prop_A = proportions) %>%
  left_join(cell_counts,
            by = c("SampleB" = "snRNA",
                   "cell_id")) %>%
  dplyr::rename(cell_prop_B = proportions) %>%
  dplyr::mutate(min_prop = ifelse(cell_prop_A <= cell_prop_B, 
                                  cell_prop_A, cell_prop_B)) %>%
  dplyr::select(cell_id, SampleA, SampleB, jsd, min_prop)


# Sample divergences
sample_divergences <- cell_type_divergences %>%
  dplyr::mutate(weighted_jsd = min_prop * jsd) %>%
  group_by(SampleA, SampleB) %>%
  summarize(patient_jsd = sum(weighted_jsd)) %>%
  pivot_wider(names_from = SampleA,
              values_from = patient_jsd)

# Multidimensional scaling 
get_JSD_mds <- function(divergences) {
  rows <- divergences[[1]]
  sample_divergences <- as.matrix(divergences[, -1])
  rownames(sample_divergences) <- rows
  
  mds_fit <- as.data.frame(cmdscale(as.dist(sample_divergences), 
                                    eig=TRUE, k=2)$points) %>%
    rownames_to_column("snRNA") %>% 
    left_join(condition_dictionary)
  
  return(mds_fit)
}

# Plot of all profiles
all_jsd_plt <- ggplot(get_JSD_mds(sample_divergences), aes(x = V1, 
                                                           y = V2, 
                                                           color = zone,
                                                           label = patient)) + 
  geom_point(size = 2) +
  ggrepel::geom_text_repel() +
  theme_classic() +
  ggtitle("Weighted JSD (all profile)")

# Do it per cell-type

cell_type_divergences <- cell_type_divergences %>%
  dplyr::select(cell_id, SampleA, SampleB, jsd) %>%
  group_by(cell_id) %>%
  nest() %>%
  dplyr::mutate(cell_type_mds = map(data, function(x) {
    get_JSD_mds(pivot_wider(x, names_from = SampleA,
                            values_from = jsd))
  })) %>%
  dplyr::mutate(mds_plts = map2(cell_type_mds, cell_id, function(x, y) {
    ggplot(x, aes(x = V1, 
                  y = V2, 
                  color = zone,
                  label = patient)) + 
      geom_point(size = 1) +
      ggrepel::geom_text_repel() +
      theme_classic() +
      ggtitle(paste0("JSD ", y))
  }))

JSD_plts = plot_grid(all_jsd_plt, cell_type_divergences$mds_plts[[1]], cell_type_divergences$mds_plts[[2]],
                     cell_type_divergences$mds_plts[[3]], cell_type_divergences$mds_plts[[4]], cell_type_divergences$mds_plts[[5]],
                     cell_type_divergences$mds_plts[[6]], cell_type_divergences$mds_plts[[7]], cell_type_divergences$mds_plts[[8]],
                     cell_type_divergences$mds_plts[[9]], cell_type_divergences$mds_plts[[10]], cell_type_divergences$mds_plts[[11]],
                     nrow = 3, ncol = 4, align = "hv")

pdf(file = "results/Figures/JSD_distances.pdf", width = 17, height = 11)

plot(JSD_plts)

dev.off()









