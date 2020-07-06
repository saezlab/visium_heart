# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Reads and generates R object with slide information
#' 1) Identifies variable genes
#' 2) Clustering + clustree for a fixed (0.2 to 1 resolution)
#' 3) Stays with highest resolution
#' 
#' It must work regardless the species
#' 

#Define aes params
cols = c('Cardiomyocytes' = '#800000',
         'Cardiomyocytes 1' = '#800000',
         'Cardiomyocytes 2' = '#9A6324',
         'Cardiomyocytes 3' = '#808000',
         'Fibroblasts' = '#911eb4',
         'Fibroblasts 1 COL15A1+' = '#911eb4',
         'Fibroblasts 1' = '#911eb4',
         'Fibroblasts 2 SCARA5+' = '#e6beff',
         'Fibroblasts 2' = '#e6beff',
         'Fibroblasts 3' = '#f032e6',
         'Endothelial cells' = '#000075',
         'Endothelial cells 1' = '#000075',
         'Endothelial cells 2' = 'blue',
         'Endothelial cells 2 POSTN+' = 'blue',
         'Endothelial cells 3' = '#568198',
         'Endothelial cells 3 PLVAP+' = '#568198',
         'Endothelial cells 3 VEGFC+' = '#568198',
         'Endothelial cells 4' = '#469990',
         'Endothelial cells 4 SEMA3G+' = '#469990',
         'Macrophages' = '#e6194B',
         'Macrophages 1 CD163+' = '#e6194B',
         'Macrophages 2 CD11C+' = '#fabebe',
         'Pericytes' = '#f58231',
         'Pericytes EGFLAM+' = '#f58231',
         'T cells' = '#ffe119',
         'Lymphatic endothelial cells' = '#ffd8b1',
         'Adipocytes' = '#000000',
         'Neuronal cells' = '#42d4f4',
         'Erythrocytes' = '#999999',
         'Proliferating cells' = '#999999',
         'Damaged endothelial cells' = '#999999',
         'Vascular smooth muscle cells' = '#aaffc3')

cols = unique(cols)

#' Runs SPARK with default parameters: to do get variable genes
#' @param dat: Seurat visium object
#' @param p_value_thrsh = p_value_thrsh to identify genes
#' @param dir_path: path of folder that includes spatial info, check Seurat
FindSparkFeatures = function(dat, p_value_thrsh = p_value_thrsh, verbose = FALSE){
  require(SPARK)
  
  rawcount = as.matrix(dat@assays$Spatial@counts)
  
  info = dat@images$slice1@coordinates[colnames(rawcount),c(2,3)]
  
  colnames(info) = c("x","y")
  
  spark = CreateSPARKObject(counts=rawcount, 
                            location=info[,1:2],
                            percentage = 0.1, 
                            min_total_counts = 10)
  
  spark@lib_size = apply(spark@counts, 2, sum)
  
  spark = spark.vc(spark, 
                   covariates = NULL, 
                   lib_size = spark@lib_size, 
                   num_core = 4,
                   verbose = verbose)
  
  spark = spark.test(spark, 
                     check_positive = T, 
                     verbose = verbose)
  
  head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
  
}

#' @param dir_path: path of folder that includes spatial info, check Seurat
#' @param var_features: how to identify variable features (top 2000 Seurat or SPARK)
#' @param p_value_thrsh: SPARK p-val threshold 
#' @param out_dir: path of folder that will contain the output results
#' @param resolution: resolution for clustering
#' @return A Seurat object with clusters as identities

process_visium = function(dir_path,
                          var_features = "seurat", 
                          p_value_thrsh = NULL,
                          out_dir = "./default",
                          resolution = 1,
                          verbose = FALSE){
  require(Seurat)
  require(clustree)
  require(cowplot)

  ## Make out folder
  system(paste0("mkdir ",out_dir))
  
  ## Define out paths
  clustree_out = paste0(out_dir,"/clustree_out.pdf")
  clust_out = paste0(out_dir,"/clustering_out.pdf")
  qc_out = paste0(out_dir,"/qc_out.pdf")
  
  ## <<Gene_expression>>
  visium_slide = Load10X_Spatial(data.dir = dir_path)
  visium_slide = SCTransform(visium_slide, assay = "Spatial", verbose = verbose)
  
  ## Basic QC plots
  qc_space = Seurat::SpatialFeaturePlot(visium_slide, features = c("nCount_Spatial",
                                                                   "nFeature_Spatial"))
  
  qc_violin = Seurat::VlnPlot(visium_slide, features = c("nCount_Spatial",
                                                         "nFeature_Spatial"))
  pdf(file = qc_out, width = 14, height = 14)
  
  plot(qc_space)
  plot(qc_violin)
  
  dev.off()
  
  ## <<Identify variable genes for clustering>>
  #To do: add SPARK wrapper
  if(var_features == "seurat"){
    
    visium_slide = FindVariableFeatures(visium_slide, selection.method = "vst", 
                       nfeatures = 2000)
    
    var_feat = VariableFeatures(visium_slide)
    
  } else if(var_features == "spark"){ #For cluster, in local it takes too much time
    var_feat = FindSparkFeatures(visium_slide, p_value_thrsh = p_value_thrsh)
  }
  
  DefaultAssay(visium_slide) = "SCT"
  
  visium_slide = RunPCA(visium_slide, 
                        assay = "SCT", 
                        verbose = verbose,
                        features = var_feat)
  
  visium_slide = FindNeighbors(visium_slide, reduction = "pca", dims = 1:30)
  # Granularity for clustree
  visium_slide = FindClusters(visium_slide, resolution = seq(from=0.1, to=1, by=0.2))
  
  pdf(file = clustree_out, height = 13,width = 8)
  # clustree
  print(clustree(visium_slide))
  dev.off()
  
  # Define granularity, default = 1
  visium_slide = FindClusters(visium_slide, verbose = verbose, 
                              resolution = resolution)
  
  # Define aesthetics
  ident_cols = sample(cols,length(levels(Idents(visium_slide))))
  
  # Run UMAP
  visium_slide = RunUMAP(visium_slide, reduction = "pca", dims = 1:30)
  
  pdf(file = clust_out, height = 6,width = 10)
  
  p1 = DimPlot(visium_slide, reduction = "umap", label = FALSE) +
    scale_color_manual(values = ident_cols)
  p2 = SpatialDimPlot(visium_slide, label = TRUE, label.size = 0,
                      stroke = 0, label.box = F) +
    scale_fill_manual(values = ident_cols)
  
  plot(plot_grid(p1,p2,nrow = 1,align = "hv"))
  
  dev.off()
  
  return(visium_slide)
}








