# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we reduce HCA object to obtain nuclei and LV
library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(tidyverse)

hca_atlas <- readRDS("./ext_data/hca_seurat.rds")

lv_atlas <- subset(hca_atlas, subset =  region == "LV" &
                      source == "Nuclei" &
                      Used == "Yes" &
                      cell_type != "doublets")

rm(hca_atlas)

saveRDS(lv_atlas, file = "./ext_data/hca_seurat_lv.rds")

print("seurat done")

# as single cell experiment
colnames(lv_atlas@meta.data) <- gsub("[.]", "_", colnames(lv_atlas@meta.data))
lv_atlas_sce <- as.SingleCellExperiment(lv_atlas)
saveRDS(lv_atlas_sce, file = "./ext_data/hca_sce_lv.rds")

print("sce done")

# scanpy ready
writeH5AD(lv_atlas_sce, file = "./ext_data/hca_lv.h5ad")

print("h5ad done")
