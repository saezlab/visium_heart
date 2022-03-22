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
library(clustree)

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
  make_option(c("--batch_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "file that contains 3 columns: orig.ident, patient, batch"),
  make_option(c("--optimize"), 
              action = "store_true", 
              default = FALSE, 
              type = 'logical',
              help = "Find best clustering? May fail for large datasets"),
  make_option(c("--default_resolution"), 
              action ="store", 
              default = 0.5, 
              type = 'double',
              help = "clustering resolution if not optimized")
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

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
# since cell-types and niches should be greater than the number of batches

hvg_list <- map(integrated_data, function(x) {
  
  DefaultAssay(x) <- def_assay
  
  x <- FindVariableFeatures(x, selection.method = "vst", 
                            nfeatures = 3000)
  
  x@assays[[def_assay]]@var.features
  
}) %>% unlist()

hvg_list <- table(hvg_list) %>%
  sort(decreasing = TRUE)

gene_selection_plt <- hvg_list %>% enframe() %>% 
  group_by(value) %>% 
  mutate(value = as.numeric(value)) %>%
  summarize(ngenes = length(name)) %>% 
  ggplot(aes(x = value, y = ngenes)) + 
  geom_bar(stat = "identity")

gene_selection <- hvg_list[1:3000] %>% names()

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
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = gene_selection, 
         npcs = 30, 
         verbose = FALSE) 

original_pca_plt <- DimPlot(object = integrated_data, 
              reduction = "pca", 
              pt.size = .1, 
              group.by = "orig.ident")

# Add batch info
# must contain orig.ident 
batch_info <- read_csv(batch_file)
batch_covars <- colnames(batch_info)

tmp_meta <- integrated_data@meta.data %>%
  left_join(batch_info, by = "orig.ident")

integrated_data@meta.data <- bind_cols(integrated_data@meta.data, tmp_meta[, batch_covars[-1]])

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              batch_covars, 
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
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name = "umap_harmony") %>%
  RunUMAP(reduction = "pca", dims = 1:30,
          reduction.name = "umap_original")

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
  
  clustree_plt <- clustree(integrated_data, 
                           prefix = paste0(DefaultAssay(integrated_data), "_snn_res."))
  
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
  
  seq_res <- seq(0.2, 1.6, 0.2)
  
  integrated_data <- FindClusters(integrated_data,
                                  resolution = seq_res,
                                  verbose = F)
  
  clustree_plt <- clustree(integrated_data, 
                           prefix = paste0(DefaultAssay(integrated_data), 
                                           "_snn_res."))
  
  integrated_data <- FindClusters(integrated_data,
                                  resolution = default_resolution,
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
        reduction = "umap_harmony", 
        pt.size = .1, 
        group.by = "orig.ident")

umap_corrected_clustering <- DimPlot(object = integrated_data, 
          reduction = "umap_harmony", 
          pt.size = .1, 
          group.by = "opt_clust_integrated")

umap_sample <- DimPlot(object = integrated_data, 
                                 reduction = "umap_original", 
                                 pt.size = .1, 
                                 group.by = "orig.ident")

umap_clustering <- DimPlot(object = integrated_data, 
                                     reduction = "umap_original", 
                                     pt.size = .1, 
                                     group.by = "opt_clust_integrated")

pdf(file = out_fig_file, height = 10, width = 12)

print(gene_selection_plt)
print(original_pca_plt)
print(corrected_pca_plt)
print(umap_sample)
print(umap_corrected_sample)
print(clustree_plt)
print(umap_clustering)
print(umap_corrected_clustering)

dev.off()

# Give reductions to ease future analysis

reductions_list <-  list(meta_data = integrated_data@meta.data,
                         reduction = integrated_data@reductions[["umap_harmony"]])

saveRDS(reductions_list,
        file = gsub("[.]rds", "_umap.rds", out_file))

