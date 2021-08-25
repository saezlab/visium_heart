# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we select from a single integrated object
#' cell-types of interest. Usually a single cell is expected
#' 
#' Once a selection is done, a whole integration process
#' and optimization of clustering is done to
#' find cell states.
#' 

library(optparse)
library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds object"),
  make_option(c("--out_fig_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the pdf object"),
  make_option(c("--cell_class"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "cell classes to consider (separated by commas)"),
  make_option(c("--class_label"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "column in data_path that contains the selected classes"),
  make_option(c("--start_res"), 
              action= "store", 
              default = 0.2, 
              type = 'double',
              help = "inital resolution for optimization")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Cell types to include ----------------------------------------------------------------------------
cell_type <- unlist(strsplit(cell_class, ","))

# Read object to subset -----------------------------------------------------
integrated_data <- readRDS(data_path)
true_ix <- integrated_data@meta.data[, class_label] %in% cell_type

# subset based on filtering and quickly get characteristic profile
integrated_data <- integrated_data[ , true_ix]

# subset the object in 2 main categories: batch and patient

integrated_data <- SplitObject(integrated_data, split.by = "batch")

# First get all the most variable genes per patient per batch

hvg_batch <- map(integrated_data, function(batch_x) {
  
  batch_x <- SplitObject(batch_x, split.by = "patient")
  
  hvg_list <- map(batch_x, function(x) {
    
    DefaultAssay(x) <- "RNA"
    
    x <- FindVariableFeatures(x, selection.method = "vst", 
                              nfeatures = 3000)
    
    x@assays[["RNA"]]@var.features
    
  }) %>% unlist()
  
  hvg_list <- table(hvg_list) %>%
    sort(decreasing = TRUE)
  
  gene_selection <- hvg_list[1:3000] %>% names()
  
})

hvg_list <- unlist(hvg_batch) %>% table() %>% sort(decreasing = T)

# Genes that are stable among batches
gene_selection <- hvg_list[hvg_list == 2] %>% names()

# Re-integrate data

integrated_data <- reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = gene_selection, 
         npcs = 30, 
         verbose = FALSE) 

# Original PCs
original_pca_plt <- DimPlot(object = integrated_data, 
                            reduction = "pca",
                            pt.size = .1,
                            group.by = "orig.ident")


# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              c("orig.ident", "patient", "batch"), 
                              plot_convergence = TRUE,
                              max.iter.harmony = 20)

# Corrected dimensions -----------------------
corrected_pca_plt <- DimPlot(object = integrated_data, 
                             reduction = "harmony", 
                             pt.size = .1, 
                             group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)

# Clustering and optimization -------------------------
print("Optimizing clustering")

integrated_data <- FindNeighbors(integrated_data, 
                                 reduction = "harmony", 
                                 dims = 1:30)

seq_res <- seq(start_res, 1, 0.2)

# Delete previous clustering
integrated_data@meta.data <- integrated_data@meta.data[, !grepl("RNA_snn_res",
                                                               colnames(integrated_data@meta.data))]

# Create new clustering
integrated_data <- FindClusters(integrated_data,
                                resolution = seq_res,
                                verbose = F)

# Optimize clustering ------------------------------------------------------
cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                   method = "euclidean")

cluster_info <- integrated_data@meta.data[,grepl("RNA_snn_res",
                                                 colnames(integrated_data@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)

silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  mean(si[, 'sil_width'])
})

integrated_data[["opt_state"]] <- integrated_data[[names(which.max(silhouette_res))]]

Idents(integrated_data) = "opt_state"

# Save object ------------------------------------------------------

saveRDS(integrated_data, file = out_path)

# Print QC file ------------------------------------------------------

umap_corrected_sample <- DimPlot(object = integrated_data, 
                                 reduction = "umap", 
                                 pt.size = .1, 
                                 group.by = "orig.ident")

umap_corrected_clustering <- DimPlot(object = integrated_data, 
                                     reduction = "umap", 
                                     pt.size = .1, 
                                     group.by = "opt_state")

# Print proportions ------------------------------------------------------

state_prop <- integrated_data@meta.data %>%
  group_by(orig.ident) %>%
  mutate(n_cells = length(orig.ident)) %>%
  ungroup() %>%
  group_by(orig.ident, opt_state) %>%
  mutate(n_state = length(opt_state)) %>%
  dplyr::select(orig.ident,n_cells, opt_state, n_state) %>%
  unique()


cell_counts <- state_prop %>%
  ungroup() %>%
  dplyr::select(orig.ident, n_cells) %>%
  unique() %>%
  ggplot(aes(x = n_cells,
             y = orig.ident)) + 
  geom_bar(stat="identity") +
  theme(legend.position = "bottom") +
  ylab("") + xlab("n_cells")   

prop_overview <- ggplot(state_prop, 
                        aes(fill = opt_state, 
                            x = n_state, 
                            y = orig.ident)) + 
  geom_bar(position="fill", 
           stat="identity") +
  theme(legend.position = "bottom") +
  ylab("") + xlab("proportions")   


pdf(file = out_fig_file, height = 10, width = 12)

print(original_pca_plt)
print(corrected_pca_plt)
print(umap_corrected_sample)
print(umap_corrected_clustering)
print(cell_counts)
print(prop_overview)

dev.off()























