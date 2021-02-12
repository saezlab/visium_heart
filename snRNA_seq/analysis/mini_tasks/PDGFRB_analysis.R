# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Generate UMAPs for mouse model

library(tidyverse)
library(Seurat)

# Read data -----------------------------------------------------------------------------------
integrated_data <- readRDS("/net/data.isilon/ag-saez/bq_shared/scellMI/processed_RNA/integrated/scRNA.integrated.rds") 
naba <- readRDS("/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/markers/NABAgsets.rds") 

# Get module scores function -----------------------------------------------------------------------------------
#' Generates a module score matrix
#' 
#' @param visium_slide: Seurat object
#' @param MS_regulon: list of gene sets
#' @return a matrix with module scores ready to be added as an assay

getTF_matrix_MS <- function(visium_slide, MS_regulon) {
  
  names_vect <- gsub("[.]","_",names(MS_regulon))
  names_vect <- gsub("-","_",names_vect)
  
  tf_act_mat <- AddModuleScore(visium_slide,
                              features = MS_regulon,
                              name = paste0(names_vect,"__"))
  
  tf_act_mat <- tf_act_mat@meta.data
  
  cell_ids <- rownames(tf_act_mat)
  calculated_regulons <- colnames(tf_act_mat)[grepl("__", colnames(tf_act_mat))]
  
  tf_act_mat <- tf_act_mat[,calculated_regulons]
  
  colnames(tf_act_mat) <- unlist(map(strsplit(colnames(tf_act_mat),split = "__"), function(x) x[1]))
  
  rownames(tf_act_mat) <- cell_ids
  
  tf_act_mat <- t(as.matrix(tf_act_mat))
  
  return(tf_act_mat)
  
}

# Main ------------------------------

# Add ECM scores ----------

DefaultAssay(integrated_data) <- "RNA"

integrated_data[['ecm_scores']] <- CreateAssayObject(data = getTF_matrix_MS(visium_slide = integrated_data,
                                                                        MS_regulon = naba))

pdf(file = "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/mini_tasks/PDGFRB_analysis.pdf", 
    height = 10, width = 12)

celltype_UMAP <- DimPlot(integrated_data, group.by = "top_annotation")
plot(celltype_UMAP)

PDGFRB_UMAP <- FeaturePlot(integrated_data, features = "PDGFRB")
plot(PDGFRB_UMAP)

DefaultAssay(integrated_data) <- "ecm_scores"

for(naba_set in rownames(integrated_data)) {
  plot(FeaturePlot(integrated_data, features = naba_set))
}

dev.off()





















































