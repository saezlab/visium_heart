# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Perform differential expression analysis of all cell-types

library(scater)
library(tidyverse)
library(rdist)
source("./analysis/utils/pseudobulk_utils.R")

# Defining functions

all_mats <- readRDS("./processed_snrnaseq/integration/psxpat_integrated_rnasamples_filt.rds")

pat_anns <- readRDS("./markers/snrna_patient_anns_revisions.rds")[, c("patient_id", "patient_group")] %>%
  unique()

marker_list_cts <- readRDS("./markers/pb_ct_marker_list.rds")

get_de_genes <- function(ct, gex_mat, pat_anns) {
  
  print("Analyzing")
  print(ct)
  
  exclude_cts <- names(marker_list_cts)
  
  exclude_cts <- exclude_cts[!(exclude_cts %in% ct)]
  
  marker_list_cts <- marker_list_cts[exclude_cts] %>%
    unlist() %>% 
    unique()
  
  # First get the samples that exist in the matrix
  
  gex_meta <- pat_anns %>%
    dplyr::filter(patient_id %in% colnames(gex_mat))
  
  groups2compare <-  gex_meta %>%
    group_by(patient_group) %>%
    summarise(n_samples = n()) %>%
    arrange(patient_group) %>%
    dplyr::filter(n_samples >= 3) %>%
    pull(patient_group)
  
  contrast_def <- combn(groups2compare,2) %>%
    t() %>%
    as_tibble() %>%
    dplyr::rename("groupA" = V1, 
                  "groupB" = V2) %>%
    dplyr::mutate(de_res = map2(groupA, groupB, 
                                run_edgeR, meta = gex_meta,
                                gex = gex_mat,
                                marker_list_cts = marker_list_cts))
}


run_edgeR <- function(groupA, groupB, meta, gex, marker_list_cts) {
  
  print("Comparing")
  print(paste0(groupA, groupB))
  
  meta_data <- meta %>%
    dplyr::filter(patient_group %in% c(groupA, groupB))
  
  sel_genes <- !(rownames(gex) %in% marker_list_cts)
  
  pb_gex <- gex[sel_genes, meta_data$patient_id]
  
  # Filter expression matrix for contrast of interest
  dat <- DGEList(pb_gex, samples = DataFrame(meta_data))
  
  keep <- filterByExpr(dat, group = meta_data$patient_group)
    
  dat <- dat[keep,]
    
  dat <- calcNormFactors(dat)
    
  design <- model.matrix(~factor(patient_group), dat$samples)
  
  dat <- estimateDisp(dat, design)
    
  fit <- glmQLFit(dat, design, robust=TRUE)
    
  res <- glmQLFTest(fit, coef=ncol(design))
    
  de_res <- topTags(res, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column("gene")

}

# MAIN

all_mats <- all_mats %>%
  mutate(de_res = map2(cell_type, count_matrix, get_de_genes, pat_anns = pat_anns))

de_res <- all_mats %>%
  dplyr::select(cell_type, de_res) %>%
  unnest() %>%
  unnest()

saveRDS(de_res, "./results/sample_comparison/all_cts_de_analysis.rds")
