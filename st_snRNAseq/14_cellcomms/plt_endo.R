# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we process results from the
#' Endo -> cts interactions
#' 
source("./analysis/utils/liana_utils.R")
library(liana)

# Get results from myofibroblast and myeloid

liana_res <- readRDS("./results/cell_comms/Endo/liana_Endo.rds")

cpdb <- liana_res$cellphonedb %>%
  dplyr::select(source, target, ligand, receptor, lr.mean)

liana_res <- liana_res %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001)) %>%
  dplyr::left_join(cpdb, by = c("source", "target", "ligand", "receptor"))


# Get the top 10 ligands 
states <- c( "Arterial_Endo", "Lymphatic_Endo", "Capillary_Endo", "Venous_Endo", "Endocardial_Endo")
cts <- c("Fib", "CM", "vSMCs", "PC")

walk(states, function(state) {
  
  print(state)
  
  liana_outs(source_groups = cts,
             target_groups = state,
             top = 10,
             file_alias = paste0(state, "_", "rec", "_"),
             max_interactions = 2,
             filter_marker_genes = T,
             out_dir = "./results/cell_comms/Endo/")
  
  liana_outs(source_groups = state,
             target_groups = cts,
             top = 10,
             file_alias = paste0(state, "_", "snd", "_"),
             max_interactions = 2,
             filter_marker_genes = T,
             out_dir = "./results/cell_comms/Endo/")
  
})



