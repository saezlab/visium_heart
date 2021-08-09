# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run the basic Seurat pipeline for spatial transcriptomics
#' folder
#' |
#' sample---
#'          |
#'          ---spatial
#'          ---filtered_feature_bc_matrix.h5
#' 1) Read 10x files
#' 2) Normalization with cpm and SCT
#' 3) Funcomics
#' 4) Clustering
#' 5) QC plots

library(optparse)
library(tidyverse)
library(Seurat)
library(cluster)
library(cowplot)
source("./analysis/utils/funcomics.R")

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
  make_option(c("--force_filter"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "shall we filter SCT matrix for MT-genes and filter spots based on ncounts?")
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
  slide_files <- paste0(path,
                        sample_names,
                        "/outs")
  } else {
  sample_names <- id
  slide_files <- path
}

# Make data frame of parameters ------------------------------------------------------------------------
param_df <- tibble(sample_name = sample_names,
                   slide_file = slide_files,
                   out_file = paste0(out_path, sample_names, ".rds"),
                   out_fig_file_0 = paste0(out_fig_path, sample_names, "_noisyspots", ".pdf"),
                   out_fig_file_a = paste0(out_fig_path, sample_names, "_slideqc", ".pdf"),
                   out_fig_file_b = paste0(out_fig_path, sample_names, "_slideclustering", ".pdf"))

# Automatic processing of each data set --------------------------------

process_data_visium <- function(sample_name, 
                                slide_file, 
                                out_file, 
                                out_fig_file_0,
                                out_fig_file_a,
                                out_fig_file_b) {
  set.seed(17)
  
  print("Reading data")
  
  print(sample_name)
  
  # We load all spots even the ones that are not over tissue
  sample_seurat <- Seurat::Load10X_Spatial(data.dir = slide_file,
                                           filename = "raw_feature_bc_matrix.h5",
                                           filter.matrix = F)
  
  coldata <- GetTissueCoordinates(sample_seurat,
                                   cols = c("row", "col", "tissue"),
                                   scale = NULL)
  
  sample_seurat$tissue <- coldata[rownames(sample_seurat@meta.data), "tissue"]
  
  sample_seurat$tissue <- ifelse(sample_seurat$tissue == 1, 
                                 "on_tissue", 
                                 "not_on_tissue")
  
  # We do a comparison of profiles between spots in tissue and not in tissue
  
  rm(coldata)
  
  tissue_qc <- sample_seurat@meta.data %>%
    select(-orig.ident) %>%
    pivot_longer(-tissue, 
                 names_to = "qc_features",
                 values_to = "counts") %>%
    ggplot(aes(x = qc_features, fill = tissue, y = counts)) +
    geom_violin() +
    facet_grid(. ~ qc_features, scales = "free")
  
  # Filter useful spots 
  
  sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")
  
  # Continue with analysis
  
  sample_seurat[["orig.ident"]] <- sample_name
  
  # Get mitochondrial genes percentage ------------------------------------------------
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, 
                                                        pattern = "^MT-")
  # QC relationships -------------------------------------------------------------------
  qc_p1 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_Spatial, y = nFeature_Spatial)) +
    geom_point() +
    theme_classic() +
    ggtitle(paste0("nspots ", ncol(sample_seurat)))
  
  qc_p2 <- sample_seurat@meta.data %>%
    ggplot(aes(x = nCount_Spatial, y = percent.mt)) +
    geom_point() +
    theme_classic() +
    ggtitle(paste0("nspots ", ncol(sample_seurat)))
  
  qc_panel <- cowplot::plot_grid(qc_p1, qc_p2, ncol = 2, align = "hv")
  
  slide_qc_p <- SpatialFeaturePlot(sample_seurat,
                                   features = c("nCount_Spatial", 
                                                "nFeature_Spatial", 
                                                "percent.mt"),
                                   ncol = 3)
  
  qc_panel_a <- cowplot::plot_grid(qc_panel, slide_qc_p, 
                                 nrow = 2, ncol = 1, 
                                 rel_heights = c(0.5, 0.5))
  
  # Process the data --------------------------------------------------------------------
  
  # Exclude MT and ribosomal genes ------------------------------------------------------
  
  if(force_filter) {
    
    #Ribosomal and mitochondrial genes are taken out
    mt_genes <- row.names(sample_seurat)[grepl("^MT-", row.names(sample_seurat))]
    rps_genes <- row.names(sample_seurat)[grepl("^RPS", row.names(sample_seurat))]
    mrp_genes <- row.names(sample_seurat)[grepl("^MRP", row.names(sample_seurat))]
    rpl_genes <- row.names(sample_seurat)[grepl("^RPL", row.names(sample_seurat))]
    rb_genes <- c(rps_genes, mrp_genes, rpl_genes)
    
    sample_seurat <- sample_seurat[!rownames(sample_seurat) %in% c(rb_genes, mt_genes), ]
    
    #Genes expressed in less that 10 spots are filtered
    sample_seurat <- sample_seurat[rowSums(GetAssayData(sample_seurat, assay = "Spatial") > 0) > 10, ]
    
    # Then I recalculate the number of genes and reads per spot
    sample_seurat$nFeature_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial") > 0)
    sample_seurat$nCount_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial"))
    
    qc_p1_filt <- sample_seurat@meta.data %>%
      ggplot(aes(x = nCount_Spatial_filt, y = nFeature_Spatial_filt)) +
      geom_point() +
      theme_classic() +
      geom_vline(xintercept = 300) +
      geom_hline(yintercept = 500) +
      ggtitle(paste0("nspots ", ncol(sample_seurat)))
    
    qc_p2_filt <- sample_seurat@meta.data %>%
      ggplot(aes(x = nFeature_Spatial_filt, y = percent.mt)) +
      geom_point() +
      theme_classic() +
      geom_vline(xintercept = 300) +
      ggtitle(paste0("nspots ", ncol(sample_seurat)))
    
    #Assumes that you have at least a single cell in a spot (same as HCA)
    sample_seurat <- subset(sample_seurat, 
                            subset = nFeature_Spatial_filt > 300 & 
                              nCount_Spatial_filt > 500)
    
    qc_panel_filt <- cowplot::plot_grid(qc_p1_filt, qc_p2_filt, ncol = 2, align = "hv")
  }
  
  # SCT transform normalization ---------------------------------------------------------
  sample_seurat <- SCTransform(sample_seurat, 
                               assay = "Spatial", 
                               verbose = FALSE)
  # cpm normalization ---------------------------------------------------------
  sample_seurat <- NormalizeData(sample_seurat, 
                                 normalization.method = 'LogNormalize', 
                                 scale.factor = 10000, 
                                 verbose = FALSE)
  
  # dimensionality reduction --------------------------------------------------
  DefaultAssay(sample_seurat) <- "SCT"
  
  sample_seurat <- ScaleData(sample_seurat, 
                             verbose = FALSE, 
                             features = rownames(sample_seurat)) %>%
    RunPCA() %>%
    RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)
  
  # Optimize clustering --------------------------------------------------------------------
  print("Optimizing clustering")
  
  sample_seurat <- FindNeighbors(sample_seurat, reduction = "pca", dims = 1:30)
  
  seq_res <- seq(0.5, 1.5, 0.1)
  
  sample_seurat <- FindClusters(sample_seurat, 
                                resolution = seq_res,
                                verbose = F)
  
  # Optimize ---------------------------------------------------------------------------------

  cell_dists <- dist(sample_seurat@reductions$pca@cell.embeddings,
                     method = "euclidean")
  
  cluster_info <- sample_seurat@meta.data[,grepl(paste0(DefaultAssay(sample_seurat), "_snn_res"),
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
  
  
  optm_res <- names(which.max(silhouette_res))
  sample_seurat[["opt_clust"]] <- sample_seurat[[optm_res]]
  
  # Reduce meta-data -------------------------------------------------------------------------
  spam_cols <- grepl(paste0(DefaultAssay(sample_seurat), "_snn_res"),
                     colnames(sample_seurat@meta.data)) |
    grepl("seurat_clusters",colnames(sample_seurat@meta.data))

  sample_seurat@meta.data <- sample_seurat@meta.data[,!spam_cols]
  
  # Plot final cluster resolution --------------------------------------------------------------------
  
  final_embedding <- DimPlot(sample_seurat, group.by = "opt_clust") +
    ggtitle(paste0("n spots ", 
                   ncol(sample_seurat), 
                   " ",
                   optm_res))
  
  spatial_embedding <- SpatialDimPlot(sample_seurat,
                                      group.by = "opt_clust",
                                      label = TRUE, 
                                      label.size = 0,
                                      stroke = 0, 
                                      label.box = F)
  
  qc_panel_b <- cowplot::plot_grid(final_embedding, spatial_embedding, ncol = 2)
  qc_panel_b_bis <- FeaturePlot(sample_seurat, features = c("percent.mt",
                                                            "nCount_Spatial",
                                                            "nFeature_Spatial"),
                                ncol = 3)
  
  print("Adding funcomics")
  
  # Add funcomics assays -
  # TF acts controlled by confidence lbls
  # pathway acts (progeny) with the top param
  sample_seurat <- add_funcomics(visium_slide = sample_seurat,
                                   species = "human",
                                   confidence_lbls = c("A","B","C","D"),
                                   top = 500,
                                   marker_dictionary = NULL,
                                   module_name = NULL,
                                   verbose = FALSE,
                                   assay = "SCT")
  
  # Save data
  print("saving object")
  Idents(sample_seurat) = "opt_clust"
  saveRDS(sample_seurat, file = out_file)
  
  # Plot QC files
  print("plotting qc features")
  
  pdf(file = out_fig_file_0, width = 10, height = 7)
  plot(tissue_qc)
  dev.off()
  
  pdf(file = out_fig_file_a, width = 16, height = 17)
  plot(qc_panel_a)
  dev.off()
  
  pdf(file = out_fig_file_b, width = 16, height = 8)
  plot(qc_panel_filt)
  plot(qc_panel_b)
  plot(qc_panel_b_bis)
  dev.off()
  
  print("finished")
}

# Main -----------------
pwalk(param_df, process_data_visium)












