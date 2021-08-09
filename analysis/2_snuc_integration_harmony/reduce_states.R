# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Identify per cell-type relevant cell-states:
#' First identify cell states that are very unique to a patient and filter them out
#' The justification is that there's no need to provide those states for deconvolution
#' 
#' Then group states based on their pseudobulk correlation using the most
#' variable genes
library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)
library(edgeR)
library(ComplexHeatmap)

source("./analysis/utils/pseudobulk_utils.R")

# Our atlas, we could expect different grouping vars:
# We filter genes that are lowly expressed
mi_atlas <- readRDS("./visium_results_manuscript/integration/ps_integrated_data_wstates.rds")[[1]]
mi_atlas_mats <- mi_atlas[names(mi_atlas) != "annotations"]
mi_atlas_mats <- map(mi_atlas_mats, assay)  %>%
  enframe("ann_level") %>%
  dplyr::mutate(origin = "nuclei",
                atlas = "mi") %>%
  dplyr::mutate(filt_mats = map(value, edgeR_filtering, min.count = 5)) %>%
  dplyr::mutate(norm_mats = map(filt_mats, cpm_norm))

# Get cell_type dict ---------------------------
scell_meta <- mi_atlas[["annotations"]]

cell_dictionary <- scell_meta %>%
  dplyr::select(cell_type, deconv_col) %>%
  unique() %>%
  dplyr::mutate(cell_type = ifelse(deconv_col == "VSMCs",
                                   "VSMCs", cell_type)) %>%
  dplyr::rename(mjr_ct = cell_type,
                cell_type = deconv_col)

# Get cell_type specific matrices
separate_counts <- separate_pseudobulk(mi_atlas_mats$filt_mats[[2]],
                                       cell_state_dic = cell_dictionary)

separate_normmats <- separate_pseudobulk(mi_atlas_mats$norm_mats[[2]],
                                         cell_state_dic = cell_dictionary)

# Get cell-state proportions ----------------------
# Manually change the VSMcs
scell_props <- scell_meta %>%
  dplyr::select(orig.ident, cell_type, deconv_col) %>%
  dplyr::mutate(cell_type = ifelse(deconv_col == "VSMCs",
                                   "VSMCs", cell_type))

# First group cells by major type and filter the ones that have more than 1 state
scell_props <- scell_props %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(w_states = map(data, function(x) {
    
    n_states <- unique(x$deconv_col)
    state_status <- ifelse(length(n_states) > 1, TRUE, FALSE)
    
  })) %>%
  unnest(w_states) %>%
  dplyr::filter(w_states == TRUE)

# Calculate basic stats of the states
# The idea of this function is to filter states that
# - have low proportion in many samples
state_calculator = function(data, 
                            within_prop = 0.01,
                            atlas_prop = 0.40) {
  
  n_patients <- length(unique(data$orig.ident))
  
  # Here we calculate the proportion of each state in
  # each sample and later filter states from samples
  # that have a low proportion
  state_stats <- data %>%
    group_by(orig.ident) %>%
    mutate(n_cells_sample = length(orig.ident)) %>%
    ungroup() %>%
    group_by(orig.ident, deconv_col) %>%
    mutate(n_cells_state = length(deconv_col)) %>%
    mutate(prop_state = n_cells_state/n_cells_sample) %>%
    unique() %>%
    ungroup() %>%
    dplyr::filter(prop_state >= within_prop)
  
  # Here we calculate how many patients are per state
    
  state_stats <- state_stats %>%
    group_by(deconv_col) %>%
    summarize(prop_patients = length(orig.ident)/n_patients) %>%
    dplyr::filter(prop_patients >= atlas_prop)
  
  return(state_stats)
  
}

# Do this for every cell-type with states
scell_props <- scell_props %>%
  mutate(state_stats = map(data, state_calculator)) %>%
  unnest(state_stats)

# Now that we have this, we can correlate the pseudobulk profiles between cells

# First identify the most variable genes per cell-type

nvar_genes <- 3000

separate_normmats <- separate_normmats %>%
  dplyr::filter(major %in% scell_props$cell_type) %>%
  dplyr::mutate(hvg_mat = map(expr_mat, function(x) {
  
    cstates <- colnames(x)
    cstates <- cstates[cstates %in% scell_props$deconv_col]
    mat <- x[, cstates]
    
    var_genes <- base::apply(mat, 1, var)
    
    var_genes <- sort(var_genes, 
                      decreasing = T)
    
    return(x[names(var_genes[1:nvar_genes]), ])
    
  })) %>%
  dplyr::mutate(state_cor = map(hvg_mat, function(x) {
    
    cor_mat <- cor(x, method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column("cell_state_a") %>%
      pivot_longer(-"cell_state_a",
                   names_to = "cell_state_b") %>%
      left_join(cell_dictionary, 
                by = c("cell_state_a" = "cell_type")) %>%
      rename(mjr_a = mjr_ct) %>%
      left_join(cell_dictionary, 
                by = c("cell_state_b" = "cell_type")) %>%
      rename(mjr_b = mjr_ct) 
    
  }))

# Create a complete dictionary -------------------------------------------------
cor_df <- dplyr::select(separate_normmats,
                        state_cor) %>%
  unnest() %>%
  dplyr::select(-major)

# Define similarity between correlations
simil_cutoff <- 0.9
state_cors <- cor_df %>%
  dplyr::filter(cell_state_a != cell_state_b,
                mjr_a == mjr_b) %>%
  dplyr::filter(cell_state_a %in% scell_props$deconv_col,
                cell_state_b %in% scell_props$deconv_col) %>%
  dplyr::mutate(identical_states = ifelse(value >= simil_cutoff, 1, 0)) %>%
  group_by(mjr_a) %>%
  nest() %>%
  dplyr::mutate(simil_mat_plt = map(data, function(x) { 
    
    ggplot(x, aes(x = cell_state_a,
                  y = cell_state_b,
                  fill= identical_states)) +
      geom_tile() +
      theme(axis.text.x = element_text(angle = 90))
      
    }))
  
pdf("./visium_results_manuscript/ct_data/state_ps_correlation.pdf")

walk(state_cors$simil_mat_plt, plot)

dev.off()


saveRDS(list("patient_filter" = scell_props,
             "correlation_filter" = state_cors),
        file = "./visium_results_manuscript/ct_data/collapsing_data.rds")
