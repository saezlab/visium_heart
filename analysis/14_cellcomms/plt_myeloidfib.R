# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we process results from the
#' Myeloid/Fibroblast interactions
#' 

source("./analysis/utils/liana_utils.R")
library(liana)

# Get results from myofibroblast and myeloid

liana_res <- readRDS("./results/cell_comms/FibMyeloid/liana_Fib_Myeloid.rds")

cpdb <- liana_res$cellphonedb %>%
  dplyr::select(source, target, ligand, receptor, lr.mean)

liana_res <- liana_res %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001)) %>%
  dplyr::left_join(cpdb, by = c("source", "target", "ligand", "receptor"))

# Get the top 10 ligands 
myeloid <- c("SPP1_Macrophages", "LYVE_FOLR_Macrophages", 
             "LYVE_PLTP_Macrophages", "CCL18_Macrophages", 
             "Monocytes")

fibroblasts <- c("Fib_0", "Myofib", "Fib_SCARA5")

# Run fibroblast focused analysis

walk(fibroblasts, function(state) {
  
  print(state)
  
  liana_outs(source_groups = myeloid,
             target_groups = state,
             top = 10,
             file_alias = paste0(state, "_", "rec", "_"),
             max_interactions = 2,
             filter_marker_genes = T, 
             out_dir = "./results/cell_comms/Fib/")
  
  liana_outs(source_groups = state,
             target_groups = myeloid,
             top = 10,
             file_alias = paste0(state, "_", "snd", "_"),
             max_interactions = 2,
             filter_marker_genes = T,
             out_dir = "./results/cell_comms/Fib/")
  
})

# Run myeloid focused analysis

walk(myeloid, function(state) {
  
  print(state)
  
  liana_outs(source_groups = fibroblasts,
             target_groups = state,
             top = 10,
             file_alias = paste0(state, "_", "rec", "_"),
             max_interactions = 3,
             filter_marker_genes = T,
             out_dir = "./results/cell_comms/Myeloid/")
  
  liana_outs(source_groups = state,
             target_groups = fibroblasts,
             top = 10,
             file_alias = paste0(state, "_", "snd", "_"),
             max_interactions = 3,
             filter_marker_gene = T,
             out_dir = "./results/cell_comms/Myeloid/")
  
})

