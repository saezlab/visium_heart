# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlation of pseudobulk profiles between atlases

library(Seurat)
library(harmony)
library(scater)
library(philentropy)
library(ComplexHeatmap)
source("./analysis/utils/funcomics.R")
source("./analysis/utils/dea.R")

ct_data <- readRDS("./visium_results_manuscript/ct_data/cardiomyocytes_states.rds")

ct_data <- getTF_matrix_MS(visium_slide = ct_data,
                           MS_regulon = hallmarks,
                           assay = "RNA",
                           module_name = "hallmarks")

deg <- readRDS("./visium_results_manuscript/ct_data/fibroblasts/fibroblasts_dea.rds")
deg[["hallmarks"]] <- find_allfeat(visium_slide = ct_data, 
                                  assays_collection = "hallmarks",
                                  logfc.threshold = 0.05,
                                  only.pos = T)[[1]]
deg <- deg %>%
  enframe() %>%
  mutate(duplicated_vars = map(value, function(x) { 
    
    x$gene[duplicated(x$gene)]
    
  })) %>%
  mutate(value = map2(value, duplicated_vars, function(x, y) {
    dplyr::filter(x, ! gene %in% y,
                  avg_logFC > 0.25,
                  pct.1 > 0.5)
  }))

dea_hmarks <- find_allfeat(visium_slide = ct_data, 
                        assays_collection = "progeny",
                        logfc.threshold = 0.05,
                        only.pos = T)


x <- "opt_state"
pseudobulk_data <- sumCountsAcrossCells(x = as.matrix(ct_data@assays$RNA@counts),
                                        ids = ct_data@meta.data[, x])
pseudobulk_data <- assay(pseudobulk_data)
dist_mat <- philentropy::JSD(t(pseudobulk_data), est.prob = "empirical")
rownames(dist_mat) = colnames(dist_mat) <- colnames(pseudobulk_data)
dist_mat <- dist_mat ** (1/2)

corrplot(dist_mat, is.corr = FALSE, method = "color", 
         tl.cex = 0.6, tl.col = "black",order = "hclust")

Heatmap(dist_mat,
        name = "JSD distance",)

colnames(pseudobulk_data) <- paste0("state_",
                                    colnames(pseudobulk_data))








duplicated()


hallmarks <- readRDS("./markers/Genesets_Dec19.rds")
hallmarks <- hallmarks$MSIGDB_HMARKS



ct_data <- getTF_matrix_MS(visium_slide = ct_data,
                           MS_regulon = hallmarks,
                           assay = "RNA",
                           module_name = "hallmarks")

DefaultAssay(ct_data) <- "hallmarks"

ct_data <- ct_data %>%
  ScaleData(verbose = FALSE) %>% 
  RunPCA(assay = "hallmarks", 
         npcs = 30, 
         verbose = FALSE,
         features = rownames(ct_data)) 

ct_data <- RunHarmony(ct_data, 
                      "orig.ident", 
                      plot_convergence = TRUE)

ct_data <- ct_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)

ct_data <- FindNeighbors(ct_data, 
                         reduction = "pca",
                         dims = 1:30)

ct_data <- FindClusters(ct_data,
                                resolution = 0.2,
                                verbose = F)

umap_corrected_sample <- DimPlot(object = ct_data, 
                                 reduction = "umap", 
                                 pt.size = .1, 
                                 group.by = "orig.ident")

umap_corrected_clustering <- DimPlot(object = ct_data, 
                                     reduction = "umap", 
                                     pt.size = .1, 
                                     group.by = "seurat_clusters")

dea_res <- find_allfeat(visium_slide = ct_data, 
                        assays_collection = "progeny",
                        logfc.threshold = 0.05,
                        only.pos = T)


dea_res$hallmarks %>% arrange(cluster, -avg_log2FC)


FeaturePlot(ct_data,
            features = c("percent.mt", 
                         "HALLMARK-INTERFERON-ALPHA-RESPONSE"),
            reduction = "umap",)











