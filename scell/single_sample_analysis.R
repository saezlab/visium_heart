# MIT License

# Copyright (c) [2019] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate QC metrics for all studies, perform
#' filtering of cells, clustering, neighbors, clustree, UMAP,
#' generate hmap of expected marker genes and run genesorter
#' 
#' QC stats returned
#' 1) Number of cells
#' 2) Median mitochondrial percentage
#' 3) Median nFeature
#' 4) Number of cells after filtering
#' 5) Median mitochondrial percentage after filtering
#' 6) Median nFeature after filtering
#' 
#' Tables returned
#' Genesorter results
#' 
#' Figures
#' QC plots
#' NumDim
#' Clustree
#' Hmaps
#' 
#' Objects
#' Seurat object of each sample
#' 

library(Seurat)
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
library(clustree)
library(genesorteR)

source("./sc_source.R")

# Load the PBMC dataset
samples = c("CK158","CK159","CK160","CK161",
            "CK162","CK163","CK164","CK165")

samples = setNames(samples,samples)

qc_stats = lapply(samples, function(unique_sample){
  
  print(unique_sample)
  
  #Make results directory
  terminal_line = sprintf("mkdir ./results/single_sample/%s",unique_sample)
  system(terminal_line)
  
  INPUT = sprintf("%s_align_res/outs/filtered_feature_bc_matrix",
                  unique_sample)
  
  input_data = Read10X(data.dir = INPUT)
  
  # Initialize the Seurat object with the raw (non-normalized data).
  sample_seurat = CreateSeuratObject(counts = input_data, 
                                     project = unique_sample, 
                                     min.cells = 3, 
                                     min.features = 200)
  
  sample_seurat[["percent.mt"]] = PercentageFeatureSet(sample_seurat, pattern = "^MT-")
  
  pdf(file = sprintf("./results/single_sample/%s/QC_metrics_flexible.pdf",
                     unique_sample), onefile=TRUE,width = 11)
  
  print(VlnPlot(sample_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3))
  
  plot1 = FeatureScatter(sample_seurat, 
                          feature1 = "nCount_RNA", 
                          feature2 = "percent.mt")
  
  plot2 = FeatureScatter(sample_seurat, 
                          feature1 = "nCount_RNA", 
                          feature2 = "nFeature_RNA")
  
  print(CombinePlots(plots = list(plot1, plot2)))
  
  sample_meta = sample_seurat@meta.data
  
  before = summarise(sample_meta, ncells = length(orig.ident),
                     med_nCount_RNA = median(nCount_RNA),
                     med_nFeature_RNA = median(nFeature_RNA),
                     med_percent.mt = median(percent.mt),
                     time = "before") 
  
  # Filtering
  nfeature_filter = quantile(sample_meta$nFeature_RNA,1-0.005)
  mt_filter = 5 #Flexible
  
  sample_seurat = sample_seurat[, (sample_seurat$nFeature_RNA > 200) & 
                                  (sample_seurat$nFeature_RNA < nfeature_filter) & 
                                  (sample_seurat$percent.mt < mt_filter)]
  
  sample_meta = sample_seurat@meta.data
  
  after = dplyr::summarise(sample_meta, ncells = length(orig.ident),
                           med_nCount_RNA = median(nCount_RNA),
                           med_nFeature_RNA = median(nFeature_RNA),
                           med_percent.mt = median(percent.mt),
                           time = "after")
  
  # Normalization
  sample_seurat = NormalizeData(sample_seurat)
  sample_seurat = FindVariableFeatures(sample_seurat, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
  
  # Scaling - before PCA
  all_genes = rownames(sample_seurat)
  sample_seurat = ScaleData(sample_seurat, 
                            features = all_genes)
  # PCA
  sample_seurat = RunPCA(sample_seurat, 
                         features = VariableFeatures(object = sample_seurat))
  
  npcs = get_npcs(seurat_object = sample_seurat,
                  show_plot = T)
  
  # Plotting and n of PCs selection
  sample_seurat = FindNeighbors(sample_seurat, dims = 1:npcs)
  
  # Granularity for clustree
  sample_seurat = FindClusters(sample_seurat, resolution = seq(from=0.1, to=1, by=0.1))
  # clustree
  print(clustree(sample_seurat))
  
  # Final resolution
  sample_seurat = FindClusters(sample_seurat, resolution = 0.5)
  
  # UMAP
  sample_seurat = RunUMAP(sample_seurat, dims = 1:npcs)
  print(DimPlot(sample_seurat, reduction = "umap"))
  
  saveRDS(sample_seurat, file = sprintf("./results/single_sample/%s/%s_res0.5.rds",
                                        unique_sample,unique_sample))
  
  #Here marker genes
  markers_suggested = read.table("./markers/Kuppe_def.txt",
                                 sep = "\t",stringsAsFactors = F,
                                 header = T)
  
  print(DoHeatmap2(SeuratObject = sample_seurat,GSC = markers_suggested,res = NULL))
  
  dev.off()
  
  gs = sortGenes(sample_seurat@assays$RNA@scale.data,
                 Idents(sample_seurat))
  
  write.files(gs,prefix = sprintf("./results/single_sample/%s/%s",
                                  unique_sample,unique_sample))
  
  return(bind_rows(before,after))
  
})

qc_stats = bind_rows(qc_stats,.id = "sample")

write.table(qc_stats, file = "./results/single_sample/gex_qc_summary.txt",
            sep = "\t",col.names = T, row.names = F, quote = F)

