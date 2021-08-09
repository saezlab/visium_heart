# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run the basic Seurat pipeline
#' 
#' 1) Create folder of results per condition
#' 2) Create folder for individual sample
#' 3) Read 10x files
#' 4) Automatically filter double cells and mitochondrial percentage
#' 5) Automatically define the dimensions of PCA
#' 6) Different cluster granularity
#' 7) UMAP

library(optparse)
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(cluster)
library(cowplot)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--folder"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%sample(/outs)/filtered_feature_bc_matrix"),
  make_option(c("--id"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "name of the sample if single path provided"),
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "direct filtered_feature_bc_matrix file"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds objects"),
  make_option(c("--out_fig_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the fig objects"),
  make_option(c("--outs_structure"), 
              action = "store", 
              default = TRUE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%sample/outs/filtered_feature_bc_matrix")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# This is an option to give a folder and parse all the files to process --------------------------------
if(folder) {
  sample_names <- list.files(path)
  if(outs_structure) {
    slide_files <- paste0(path,sample_names,"/outs/filtered_feature_bc_matrix")
  } else {
    slide_files <- paste0(path,sample_names,"/filtered_feature_bc_matrix")
  }
} else {
  sample_names <- id
  slide_files <- path
}

# Make data frame of parameters ------------------------------------------------------------------------

param_df <- tibble(sample_name = sample_names,
                   slide_file = slide_files,
                   out_file = paste0(out_path, sample_names, ".rds"),
                   out_fig_file = paste0(out_fig_path, sample_names, ".pdf"))

# Automatic processing of each data set --------------------------------

process_data <- function(sample_name, slide_file, out_file, out_fig_file) {
  set.seed(17)
  
  print("Reading data")
  
  print(sample_name)
  
  input_data <- Read10X(data.dir = slide_file)
  
  # Initialize the Seurat object
  sample_seurat <- Seurat::CreateSeuratObject(counts = input_data, 
                                              project = sample_name, 
                                              min.cells = 10, 
                                              min.features = 300)
  rm(input_data)
  
  # Get mitochondrial genes
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, pattern = "^MT-")
  
  print("Filtering cells by gene expression and mitochondrial gene expression")
  
  # Get a less stringent quantile check (since we are taking doublets anyway)
  # I will take the 1%
  nfeature_filter = quantile(sample_seurat$nFeature_RNA,
                             1-0.01)
  
  
  # Here we make a collection of plots needed ------------------------------
  # N features and ncount
  
  filt_p1 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 300) + 
    geom_hline(yintercept = nfeature_filter) +
    geom_vline(xintercept = 500) +
    ggtitle(paste0("ncells ", ncol(sample_seurat)))
  
  filt_p2 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_RNA, y = percent.mt)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 5) + 
    geom_vline(xintercept = 500) +
    ggtitle(paste0("ncells ", ncol(sample_seurat)))
  
  # Get mitochondrial genes
  sample_seurat <- subset(sample_seurat, 
                          subset = nFeature_RNA > 200 & 
                            nFeature_RNA < nfeature_filter & 
                            percent.mt < 5 & 
                            nCount_RNA > 300)
  
  print("Identifying doublets")
  
  # Identify doublets
  doublets <- scDblFinder(sce = as.matrix(sample_seurat@assays$RNA@counts))
  sample_seurat$doublet_score <- doublets$scDblFinder.score
  sample_seurat$doublet <- doublets$scDblFinder.class
  
  print("Getting QC info")
  
  # Process the data --------------------------------------------------------------------
  sample_seurat <- NormalizeData(sample_seurat, 
                                 normalization.method = 'LogNormalize', 
                                 scale.factor = 10000, 
                                 verbose = FALSE)
  
  sample_seurat <- FindVariableFeatures(sample_seurat, 
                                        selection.method = 'vst', 
                                        nfeatures = 3000, 
                                        verbose = FALSE)
  sample_seurat <- ScaleData(sample_seurat, 
                             verbose = FALSE, 
                             features = rownames(sample_seurat))
  
  sample_seurat <- RunPCA(sample_seurat, 
                          verbose = FALSE)
  
  sample_seurat <- RunUMAP(sample_seurat, reduction = 'pca', dims = 1:30, verbose = FALSE)
  
  quality_plt <- FeaturePlot(sample_seurat, features = c("doublet_score", 
                                          "percent.mt", 
                                          "nCount_RNA",
                                          "nFeature_RNA"))
  
  quality_plt_bis <- DimPlot(sample_seurat, group.by = "doublet")
  
  # Filter and do it again --------------------------------------------------------------------
  print("Getting variable features, scaling and low dimension representations")
  
  sample_seurat <- subset(sample_seurat, 
                          subset = doublet == "singlet")
  
  sample_seurat <- FindVariableFeatures(sample_seurat, 
                                        selection.method = 'vst', 
                                        nfeatures = 3000, 
                                        verbose = FALSE)
  
  sample_seurat <- ScaleData(sample_seurat, 
                             verbose = FALSE, 
                             features = rownames(sample_seurat))
  
  sample_seurat <- RunPCA(sample_seurat, 
                          verbose = FALSE)
  
  sample_seurat <- RunUMAP(sample_seurat, 
                           reduction = 'pca', 
                           dims = 1:30, 
                           verbose = FALSE)
  
  # Generate different clustering resolutions --------------------------------------------------------------------
  print("Optimizing clustering")
  
  sample_seurat <- FindNeighbors(sample_seurat, reduction = "pca", dims = 1:30)
  
  seq_res <- seq(0.5, 1.5, 0.1)
  
  sample_seurat <- FindClusters(sample_seurat, 
                              resolution = seq_res,
                               verbose = F)
  
  # Optimize ---------------------------------------------------------------------------------
  
  cell_dists <- dist(sample_seurat@reductions$pca@cell.embeddings,
                     method = "euclidean")
  
  cluster_info <- sample_seurat@meta.data[,grepl("RNA_snn_res",
                                                 colnames(sample_seurat@meta.data))] %>%
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
  
  sample_seurat[["opt_clust"]] <- sample_seurat[[names(which.max(silhouette_res))]]
  
  # Plot final resolution
  
  final_embedding <- DimPlot(sample_seurat, group.by = "opt_clust") +
    ggtitle(paste0("n cells ", ncol(sample_seurat)))
  
  print("Generating outputs")

  # Save data
  Idents(sample_seurat) = "opt_clust"
  saveRDS(sample_seurat, file = out_file)
  
  # Plot QC files
  
  pdf(file = out_fig_file, width = 15, height = 15)
  
  plot(cowplot::plot_grid(nrow = 1, ncol = 2, filt_p1, filt_p2))
  plot(quality_plt)
  plot(quality_plt_bis)
  plot(final_embedding)
  
  dev.off()

}

# Main -----------------

pwalk(param_df, process_data)

