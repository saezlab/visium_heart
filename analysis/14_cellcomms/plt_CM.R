# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we process results from the
#' CM -> state interactions
#' 
source("./analysis/utils/liana_utils.R")
library(liana)

# Get results from myofibroblast and myeloid

liana_res <- readRDS("./results/cell_comms/CM/liana_CM.rds")

cpdb <- liana_res$cellphonedb %>%
  dplyr::select(source, target, ligand, receptor, lr.mean)

liana_res <- liana_res[which(names(liana_res) != "sca")] %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001)) %>%
  dplyr::left_join(cpdb, by = c("source", "target", "ligand", "receptor"))

# Get the top 10 ligands 
cm <- c("damaged_CM")
cts <- c("Fib", "Adipo", "Endo", "vSMCs", "Myeloid")

liana_outs(source_groups = cts,
           target_groups = cm,
           top = 10,
           file_alias = paste0(cm, "_", "rec", "_"),
           max_interactions = 2,
           out_dir = "./results/cell_comms/CM/")

liana_outs(source_groups = cm,
           target_groups = cts,
           top = 10,
           file_alias = paste0(cm, "_", "snd", "_"),
           max_interactions = 2,
           out_dir = "./results/cell_comms/CM/")


