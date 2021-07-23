# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we automatize the integration of Seurat objects using Harmony
#' 
#' It requires a folder path with processed Seurat objects from run_singleprocessing.R

library(optparse)
library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "location where samples are located"),
  make_option(c("--out_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds objects"),
  make_option(c("--out_fig_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the plots objects"),
  make_option(c("--def_assay"), 
              action= "store", 
              default = "RNA", 
              type = 'character',
              help = "default assays to integrate"),
  make_option(c("--def_assay"), 
              action= "store", 
              default = "RNA", 
              type = 'character',
              help = "default assays to integrate"),
  make_option(c("--optimize"), 
              action = "store", 
              default = TRUE, 
              type = 'logical',
              help = "Find best clustering? May fail for large datasets")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

slide_files <- list.files(path)
# Because of incompatibility with objects objects should be appended manually
slide_files_path <- set_names(paste0(path,slide_files), gsub(pattern = "[.]rds",
                                                             replacement = "",
                                                             slide_files))

integrated_data <- map(slide_files_path, readRDS)
print("You managed to load everything")
print("Object size")
print(object.size(integrated_data))

# Create merged object ---------------------------------
integrated_data <- reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

print("You managed to merge everything")
print("Object size")
print(object.size(integrated_data))

# Default assay ---------------------------------------
DefaultAssay(integrated_data) <- def_assay

# Process it before integration -----------------------
integrated_data <- integrated_data %>%
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = integrated_data@var.genes, 
         npcs = 30, 
         verbose = FALSE) 

original_pca_plt <- DimPlot(object = integrated_data, 
              reduction = "pca", 
              pt.size = .1, 
              group.by = "orig.ident")


# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              "orig.ident", 
                              plot_convergence = TRUE,
                              assay.use = def_assay,
                              max.iter.harmony = 20)

# Corrected dimensions -----------------------
corrected_pca_plt <- DimPlot(object = integrated_data, 
              reduction = "harmony", 
              pt.size = .1, 
              group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)

integrated_data <- FindNeighbors(integrated_data, 
                                 reduction = "harmony", 
                                 dims = 1:30)

if(optimize) { 
  
  # Clustering and optimization -------------------------
  print("Optimizing clustering")
  
  seq_res <- seq(0.5, 1.5, 0.1)
  
  integrated_data <- FindClusters(integrated_data,
                                  resolution = seq_res,
                                  verbose = F)
  
  # Optimize clustering ------------------------------------------------------
  cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                     method = "euclidean")
  
  
  cluster_info <- integrated_data@meta.data[,grepl(paste0(DefaultAssay(integrated_data),"_snn_res"),
                                                   colnames(integrated_data@meta.data))] %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::mutate_all(as.numeric)
  
  silhouette_res <- apply(cluster_info, 2, function(x){
    si <- silhouette(x, cell_dists)
    if(!is.na(si)) {
      mean(si[, 'sil_width'])
    } else {
      NA
    }
  })
  
  integrated_data[["opt_clust_integrated"]] <- integrated_data[[names(which.max(silhouette_res))]]
  
  Idents(integrated_data) = "opt_clust_integrated"
  
  # Reduce meta-data -------------------------------------------------------------------------
  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))
  
  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]
  
} else {
  
  print("Not Optimizing clustering")
  
  integrated_data <- FindClusters(integrated_data,
                                  resolution = 1,
                                  verbose = F)
  
  integrated_data[["opt_clust_integrated"]] <- integrated_data[["seurat_clusters"]]
  
  spam_cols <- grepl(paste0(DefaultAssay(integrated_data), "_snn_res"),
                     colnames(integrated_data@meta.data)) |
    grepl("seurat_clusters",colnames(integrated_data@meta.data))
  
  integrated_data@meta.data <- integrated_data@meta.data[,!spam_cols]
  
}

# Save object ------------------------------------------------------
saveRDS(integrated_data, file = out_file)

# Print QC file ------------------------------------------------------

umap_corrected_sample <- DimPlot(object = integrated_data, 
        reduction = "umap", 
        pt.size = .1, 
        group.by = "orig.ident")

umap_corrected_clustering <- DimPlot(object = integrated_data, 
          reduction = "umap", 
          pt.size = .1, 
          group.by = "opt_clust_integrated")

pdf(file = out_fig_file, height = 10, width = 12)

print(original_pca_plt)
print(corrected_pca_plt)
print(umap_corrected_sample)
print(umap_corrected_clustering)

dev.off()
