# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create hmaps per patient, per cell-type of interest

library(scater)
library(ComplexHeatmap)
library(tidyverse)
source("./analysis/utils/pseudobulk_utils.R")

gsets <- readRDS("./markers/gene_list_12_04_21.rds")

pseudobulk_data <- readRDS("./visium_results_manuscript/integration/ps_patients_ct.rds")

condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(snRNA, New.Ids) %>%
  dplyr::rename(orig.ident = snRNA,
                patient = New.Ids)

# First: Extract data of relevant cell-types
pseudobulk_data <- pseudobulk_data[[1]]$gex

col_meta <- colData(pseudobulk_data) %>%
  as.data.frame() %>%
  left_join(condition_dictionary)

cts <- set_names(c("fibroblasts", 
                   "cardiomyocytes", 
                   "endothelial_cells"))

ct_ps <- map(cts, function(ct) {
  
  # First, get indexes
  
  ct_ix <- which(col_meta %>% pull(cell_type) == ct)
  
  ct_dat <- assay(pseudobulk_data)[, ct_ix]
  colnames(ct_dat) <- (col_meta %>% pull(patient))[ct_ix]
  
  ct_dat <- cpm_norm(edgeR_filtering(ct_dat))
  
  ct_dat <- t(scale(t(ct_dat)))
  
})

hmaps_plts <- map(names(gsets), function(gs_name){
  
  print(gs_name)
  
  gs <- gsets[[gs_name]]
  
  map(names(ct_ps), function(ct_ps_name) { 
    
    print(ct_ps_name)
    
    mat <- ct_ps[[ct_ps_name]]
      
    genes <- gs[gs %in% rownames(mat)]
    
    hmap_plt <- ComplexHeatmap::Heatmap(mat[genes,],
                                        name = "expr",
                                        show_row_dend = F)
    
    write.table(mat[genes,], file = paste0("./visium_results_manuscript/mini_tasks/",
                                           gs_name, "_", ct_ps_name, ".txt") ,
                sep = "\t",row.names = T, col.names = T, quote = F)
    
  })
  
})

pdf("./visium_results_manuscript/mini_tasks/marker_hmaps.pdf", width = 7, height = 10)

walk(set_names(names(hmaps_plts)), function(markers_id) { 
  
  marker_hmaps <- hmaps_plts[[markers_id]]
  
  walk(set_names(names(marker_hmaps)), function(cell_id) { 
    
    draw(marker_hmaps[[cell_id]], 
         column_title = paste0(markers_id, "_", cell_id))
    
    
    })
  
  
  })

dev.off()

















