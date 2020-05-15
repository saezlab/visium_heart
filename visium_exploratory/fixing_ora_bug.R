# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here, we will analyze single slides.
#' We will map single cell markers whenever is possible
#' Data is already normalized, we need to estimate variable features
#' Create clustering and enrich cluster with gene markers and functions

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(cowplot)
library(clustree)
library(xlsx)
library(genesorteR)

source("./visium_exploratory/slide_processing.R")
gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")


# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

# Final tables per-slide

for(slide in slides_ids){
  print(slide)
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)

  dea_res = readRDS(dea_file)
  
  gene_res = dea_res$SCT %>% 
    filter(avg_logFC > 0,
           p_val_adj < 0.005) %>% 
    arrange(p_val_adj, avg_logFC) %>%
    group_by(cluster) %>%
    slice(1:100) 
  
  clusters_test = unique(gene_res$cluster)
  names(clusters_test) = clusters_test
  
  ora_results = lapply(clusters_test, function(x){
    
    gset = gene_res %>% 
      dplyr::filter(cluster == x) %>%
      dplyr::select(gene) %>% pull()
    
    GSE_analysis_res = GSE_analysis(geneList = gset,
                                    Annotation_DB = gene_sets$MSIGDB_CANONICAL)
    
    GSE_analysis_res$Cluster_ID = paste0("cluster",as.character(x),collapse = "_")
    
    return(GSE_analysis_res)
    
  }) %>% enframe("iteration") %>% unnest() %>%
    dplyr::filter(corr_p_value<0.1)
  
  dea_res[["ora_canonical"]] = ora_results
  
  saveRDS(dea_res, file = dea_file)
}



