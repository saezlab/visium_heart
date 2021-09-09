# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we perform iterative clustering of the ILR data in a Seurat object
library(compositions)
library(Seurat)
library(tidyverse)


ilr_data <- readRDS("./processed_visium/integration/integrated_slides_ct_ILR.rds")

seurat_ilr <- CreateSeuratObject(t( ilr_data$ilr_counts %>% as.data.frame() %>% as.matrix() ), 
                                 project = "SeuratProject", 
                                 assay = "ILR",
                                 in.cells = 0, min.features = 0, names.field = 1,
                                 names.delim = "_", meta.data = ilr_data$meta_data %>% column_to_rownames("row_id"))

neighbor_graph <- FindNeighbors(GetAssay(seurat_ilr, "ILR"))
seurat_ilr@graphs$ILR_snn <- neighbor_graph$snn
seurat_ilr@graphs$ILR_nn <- neighbor_graph$nn

print("Created network")

seurat_ilr <- FindClusters(seurat_ilr, resolution = seq(0.4, 1.6, 0.2))

print("Found clusters")

saveRDS(seurat_ilr, file = "./processed_visium/integration/integrated_slides_ct_ILR_srt.rds")




