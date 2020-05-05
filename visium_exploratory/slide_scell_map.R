# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here, we will analyze single slides.
#' We will map single cell markers whenever is possible
#' Data is already normalized, we need to estimate variable features
#' Create clustering and enrich cluster with gene markers and functions

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(cowplot)
library(clustree)
library(xlsx)
library(genesorteR)

source("./visium_exploratory/slide_processing.R")
gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")


# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

# Re-clustering + Wilcoxon
for(slide in slides_ids){
  
  print(slide)
  # File definition
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)
  
  clustree_file = sprintf("./results/single_slide/%s/%s_clustree.pdf",
                       slide,slide)
  
  clust_file = sprintf("./results/single_slide/%s/%s_clustering.pdf",
                          slide,slide)
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  # Slide clustering
  
  visium_slide = readRDS(slide_file)
  
  DefaultAssay(visium_slide) = "SCT"
  
  visium_slide = RunPCA(visium_slide, assay = "SCT", verbose = FALSE)
  visium_slide = FindNeighbors(visium_slide, reduction = "pca", dims = 1:30)
  # Granularity for clustree
  visium_slide = FindClusters(visium_slide, resolution = seq(from=0.1, to=1, by=0.2))
  
  pdf(file = clustree_file, height = 10,width = 8)
  # clustree
  print(clustree(visium_slide))
  dev.off()
  
  visium_slide = FindClusters(visium_slide, verbose = FALSE,resolution = 1)
  visium_slide = RunUMAP(visium_slide, reduction = "pca", dims = 1:30)
  
  
  pdf(file = clust_file, height = 4,width = 8)
  
  p1 = DimPlot(visium_slide, reduction = "umap", label = TRUE)
  p2 = SpatialDimPlot(visium_slide, label = TRUE, label.size = 3)
  
  plot(plot_grid(p1,p2,nrow = 1,align = "hv"))
  
  dev.off()
  
  # Differential feature analysis
  
  possible_assays = c("SCT","dorothea","progeny",
                      "ctscores","ECM","ctscores_match")
  
  names(possible_assays) = possible_assays
  
  possible_assays = possible_assays[possible_assays %in%
                                      Assays(visium_slide)]
  
  diff_features = lapply(possible_assays, function(x){
    
    
    DefaultAssay(visium_slide) = x
    
    vis_wilcox_markers = FindAllMarkers(visium_slide, 
                                        logfc.threshold = 0.05)
    
  })
  
  # Gene_sorter
  
  gs = sortGenes(visium_slide@assays$SCT@scale.data,
                 Idents(visium_slide))
  
  diff_features[["gs_condGeneProb"]] = data.frame(gs$condGeneProb) %>% 
    rownames_to_column("gene")
  
  diff_features[["gs_postClustProb"]] = data.frame(gs$postClustProb) %>% 
    rownames_to_column("gene")
  
  diff_features[["gs_specScore"]] = data.frame(gs$specScore) %>% 
    rownames_to_column("gene")
  
  
  saveRDS(diff_features, file = dea_file)
  
  saveRDS(visium_slide, file = slide_file)
  
  
}


# Final tables per-slide

for(slide in slides_ids){
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  final_file = sprintf("./results/single_slide/%s/%s_differential_features.xlsx",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  assays_res = names(dea_res)
  names(assays_res) = assays_res
  
  possible_assays = c("SCT","dorothea" ,"progeny",
                      "ECM","ctscores",
                      "ctscores_match")
  
  names(possible_assays) = possible_assays
  possible_assays = possible_assays[possible_assays %in%
                                      names(dea_res)]
  
  sign_tables = lapply(possible_assays, function(x){
    
    sum_table = dea_res[[x]] %>% 
      arrange(cluster,p_val_adj, avg_logFC) %>%
      dplyr::filter(avg_logFC > 0,
                    p_val_adj<=0.005)
    
    if(nrow(sum_table)>0){
    
    write.xlsx(sum_table, file = final_file, 
               sheetName = x, append = TRUE)
    }
    
  })
  
  gene_res = dea_res$SCT %>% 
    filter(avg_logFC > 0,
           p_val_adj < 0.005) %>% 
    arrange(p_val_adj, avg_logFC) %>%
    group_by(cluster) %>%
    slice(1:100) 
  
  clusters_test = unique(gene_res$cluster)
  names(clusters_test) = clusters_test
  
  ora_results = lapply(clusters_test, function(x){
    
    gset = gene_res %>% 
      dplyr::filter(cluster == x) %>%
      dplyr::select(gene) %>% pull()
    
    GSE_analysis_res = GSE_analysis(geneList = gset,
                 Annotation_DB = gene_sets$MSIGDB_CANONICAL)
    
    GSE_analysis_res$Cluster_ID = paste0("cluster",as.character(x),collapse = "_")
    
    return(GSE_analysis_res)
    
  }) %>% enframe("iteration") %>% unnest() %>%
    dplyr::filter(corr_p_value<0.1)
  
  dea_res[["ora_canonical"]] = ora_results
  
  write.xlsx(ora_results, file = final_file, 
             sheetName = "ora_sct", append = TRUE)
  
  # Gene sorter
  
  write.xlsx(dea_res[["gs_condGeneProb"]], file = final_file, 
             sheetName = "gs_condGeneProb", append = TRUE)
  
  write.xlsx(dea_res[["gs_postClustProb"]], file = final_file, 
             sheetName = "gs_postClustProb", append = TRUE)
  
  write.xlsx(dea_res[["gs_specScore"]], file = final_file, 
             sheetName = "gs_specScore", append = TRUE)

  saveRDS(dea_res, file = dea_file)
}


# Plots

for(slide in slides_ids){
  print(slide)
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  plot_file = sprintf("./results/single_slide/%s/%s_de_features.pdf",
                      slide,slide)
  
  DefaultAssay(visium_slide) = "SCT"
  
  # Dot plot of top 5 per cluster
  
  gene_summary = dea_res[["SCT"]] %>% 
    arrange(cluster,p_val_adj, avg_logFC) %>%
    dplyr::filter(avg_logFC > 0,
                  p_val_adj<=0.005) %>%
    group_by(cluster) %>%
    slice(1:5)
  
  sct_dots = DotPlot(visium_slide,features = unique(gene_summary$gene),
          assay = "SCT") + coord_flip() +
    theme(axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=8,angle = 90,
                                     vjust = 0.5),
          axis.title = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))
  
  # Dot plot of top 5 TFs
  gene_summary = dea_res[["dorothea"]] %>% 
    arrange(cluster,p_val_adj, avg_logFC) %>%
    dplyr::filter(avg_logFC > 0,
                  p_val_adj<=0.005) %>%
    group_by(cluster) %>%
    slice(1:5)
  
  dorothea_dots = DotPlot(visium_slide,
                          features = unique(gene_summary$gene),
          assay = "dorothea") + coord_flip() +
    theme(axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=8,angle = 90,
                                     vjust = 0.5),
          axis.title = element_blank(),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))
  
  # All from the rest
  possible_assays = c("progeny","ECM","ctscores",
                      "ctscores_match")
  names(possible_assays) = possible_assays
  possible_assays = possible_assays[possible_assays %in%
                                      Assays(visium_slide)]
  
  plot_summary_data = lapply(dea_res[possible_assays], function(x){
    x %>% arrange(cluster,p_val_adj, avg_logFC) %>%
      dplyr::filter(avg_logFC > 0,
                    p_val_adj<=0.005) %>%
      group_by(cluster) %>%
      slice(1:5)
  })
  
   plot_list = lapply(names(plot_summary_data), function(x){
    
     if(nrow(plot_summary_data[[x]])>0){
     
    feat_list =  unique(plot_summary_data[[x]]$gene)
    
    DotPlot(visium_slide,features = feat_list,
            assay = x) +
       coord_flip() + theme(axis.text.y = element_text(size=7),
                            axis.text.x = element_text(size=8,angle = 90,
                                                       vjust = 0.5),
                            axis.title = element_blank(),
                            legend.text = element_text(size=8),
                            legend.title = element_text(size=8))
     }else{
       NULL
     }
    
  })
  
  upper = plot_grid(sct_dots, dorothea_dots, align = "vh",
            nrow = 1)
  
  lower = plot_grid(plotlist=plot_list, align = "vh")
  
  panel_all = plot_grid(upper, lower,
                        nrow=1)
  
  pdf(file = plot_file,
      width = 19,height = 10)
  
  plot(panel_all)
  
  dev.off()
   
}











