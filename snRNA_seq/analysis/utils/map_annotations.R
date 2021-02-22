library(Seurat)
library(tidyverse)

vars_to_transfer <- "States"

data_A <- readRDS("./visium_results_manuscript/ct_data/cardiomyocytes_states.rds")
data_B <- readRDS("./visium_results_manuscript/ct_data/snRNA.Rds")

# To be sure we always use cell_ids and orig.ident
meta_son <- data_B@meta.data %>%
  rownames_to_column("cell_id")

meta_son <- meta_son[,c("cell_id","orig.ident",vars_to_transfer)]
rm(data_B)

meta_parent <- data_A@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(meta_son, by = c("cell_id","orig.ident"))


data_A <- AddMetaData(object = data_A,meta_parent[,vars_to_transfer], 
                      col.name = vars_to_transfer)

saveRDS(data_A, file = "./visium_results_manuscript/ct_data/cardiomyocytes_states.rds")
