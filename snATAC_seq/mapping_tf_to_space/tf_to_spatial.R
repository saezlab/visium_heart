# Copyright (c) [2020] [Ivan G. Costa / Zhijian Li]
# ivan.costa@rwth-aachen.de

# Spatial Mapping of HINT based TF Activities 
require(Seurat)
require(ggplot2)
require(matrixStats)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(cowplot)
library(scatterpie)
library(dorothea)
library(doParallel)
library(foreach)
registerDoParallel(cores=30)

atac_to_spatial <- c("CK166" = "157771",
                     "CK167" = "157785",
                     "CK168" = "157781",
                     "CK170" = "157772",
                     "CK171" = "157777",
                     "CK173" = "157782",
                     "CK174" = "157775")

atac_to_patient <- c("CK166" = "P1_RZ",
                     "CK167" = "P5_CR",
                     "CK168" = "P3_BZ",
                     "CK170" = "P4_CR",
                     "CK171" = "P2_BZ",
                     "CK173" = "P3_RZ",
                     "CK174" = "P2_IZ")

# accessing the cell composition of visium slide
# create assay for prediction score
ct_scores_lt_all = lapply(readRDS(file = "../results_slides/single_slide/integrated_meta.rds"),
                          function(x){
                              x = x %>% rownames_to_column("spot_id") %>%
                                  dplyr::mutate(predicted.id = gsub("[.]","_",
                                                                    gsub(" ","_", predicted.id)))
                              colnames(x) = gsub("[.]","_", colnames(x))
                              return(x)
                          })

# only cells with a matching group in scATAC and scRNA are used
cell_names_mapping <- c("Cardiomyocytes" = "Cardiomyocytes",
                        "Endothelial_1" = "Endothelial_cells_1",
                        "Endothelial_2" = "Endothelial_cells_2",
                        "Endothelial_3" = "Endothelial_cells_3",
                        "Endothelial_4" = "Endothelial_cells_4",
                        "Fibroblasts_1" = "Fibroblasts_1",
                        "Fibroblasts_2" = "Fibroblasts_2",
                        "Macrophages_1" = "Macrophages_1",
                        "Neuronal" = "Neuronal_cells",
                        "Pericytes" = "Pericytes",
                        "T_cell" = "T_cells",
                        "VSMCs" = "Vascular_smooth_muscle_cells",
                        "Lymphatic_Endothelial" = "Lymphatic_endothelial_cells")

values_all <- c()
for(sample in names(atac_to_spatial)){
    output_dir <- sprintf("%s_%s", sample, atac_to_patient[sample])
    
    # loading TF activity scores of corresponding samples 
    tf_activity_table <- read.delim(sprintf("../HINT/DiffFootprintsPerSample/%s/heatmap/TF.txt", sample))
    
    # filtering TFs with low number of binding sites
    # tf_activity_table <- tf_activity_table[tf_activity_table$Num>1000,]
    
    #creating a matrix only containg TF names and TF activity z-scores
    tf_activity = as.matrix(tf_activity_table[, 3:(dim(tf_activity_table)[2]-2)])
    #removing motif names 
    gene_name<-function(string)
    {
        toupper(unlist(strsplit(string, "\\."))[3])
    }
    rownames(tf_activity) <- unlist(lapply(as.character(tf_activity_table[,1]), gene_name))
    
    ct_scores <- as.data.frame(ct_scores_lt_all[[atac_to_spatial[sample]]])
    rownames(ct_scores) <- ct_scores$spot_id
    ct_scores <- ct_scores[, grep("prediction_score", colnames(ct_scores))]
    ct_scores$prediction_score_max <- NULL
    colnames(ct_scores) <- stringr::str_split_fixed(colnames(ct_scores), "_", 3)[, 3]
    
    # update cell type names of scATAC to spatial data
    colnames(tf_activity) <- stringr::str_replace_all(colnames(tf_activity),
                                                             cell_names_mapping)
    
    cell_scores <- as.matrix(ct_scores[colnames(tf_activity)])
    cell_scores <- cell_scores / rowSums(cell_scores)
    
    #Creating new Seurat Essay
    tf_activity_patches <- tf_activity %*% t(cell_scores)
    colnames(tf_activity_patches) <- rownames(ct_scores)
    
    # loading visium slides
    visium_slide <- readRDS(sprintf("../results_slides/single_slide/%s/%s.rds", 
                                    atac_to_spatial[sample],
                                    atac_to_spatial[sample]))
    
    visium_slide[['HINT_TF_Activity']] <- CreateAssayObject(counts = tf_activity_patches)
    DefaultAssay(visium_slide) <- 'HINT_TF_Activity'
    visium_slide <- ScaleData(visium_slide)
    visium_slide <- RunPCA(visium_slide, features = rownames(visium_slide), verbose = FALSE)

    visium_slide <- FindNeighbors(visium_slide, dims = 1:30, verbose = FALSE)
    visium_slide <- RunUMAP(visium_slide, dims = 1:30, umap.method = "uwot", metric = "cosine")
    
    # use the clusters gained by gene expression
    # visium_slide <- FindClusters(visium_slide, resolution = 0.15, verbose = FALSE)
    
    if(!dir.exists(output_dir)){
        dir.create(output_dir)
    }
    if(!dir.exists(sprintf("%s/hint_activity", output_dir))){
        dir.create(sprintf("%s/hint_activity", output_dir))
    }
    if(!dir.exists(sprintf("%s/dorothea", output_dir))){
        dir.create(sprintf("%s/dorothea", output_dir))
    }
    if(!dir.exists(sprintf("%s/expression", output_dir))){
        dir.create(sprintf("%s/expression", output_dir))
    }    
    
    # Define aesthetics
    Idents(visium_slide) <- "SCT_snn_res.1"
    
    p1 <- DimPlot(visium_slide, reduction = "umap", label = TRUE, pt.size = 0.5) + 
        NoLegend() +
        scale_color_brewer(palette = "Paired")
    p2 <- SpatialDimPlot(visium_slide, label = TRUE, label.size = 0,
                        stroke = 0, label.box = F) +
        scale_fill_brewer(palette = "Paired")
    pdf(sprintf("%s/Cluster_TF_activity.pdf", output_dir))
    print(p1)
    print(p2)
    dev.off()
    
    df_clusters <- as.data.frame(visium_slide$SCT_snn_res.1)
    colnames(df_clusters) <- c("clusters")
    cell_scores <- as.data.frame(cell_scores)
    cell_scores$clusters <- df_clusters$clusters
    
    cells_clusters <- melt(cell_scores)
    p1 <- ggplot(cells_clusters, aes(x=clusters, y=value, fill=variable)) + 
        geom_boxplot() +
        scale_fill_brewer(palette = "Paired")
        
    p2 <- ggplot(cells_clusters, aes(x=variable, y=value, fill=clusters)) + 
        geom_boxplot() + 
        scale_fill_brewer(palette = "Paired") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    pdf(sprintf("%s/Cells_Cluster.pdf", output_dir), height = 12, width = 12)
    print(p1)
    print(p2)
    dev.off()
    
    markers <- FindAllMarkers(visium_slide,
                              only.pos = TRUE,
                              logfc.threshold = 0.0,
                              verbose = TRUE)
    
    top_tfs <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    
    # generate heatmap for top 10 TFs of each cluster
    p <- DoHeatmap(visium_slide,
                   features = top_tfs$gene,
                   group.by = "SCT_snn_res.1",
                   slot = "scale.data")
    pdf(sprintf("%s/Heatmap.pdf", output_dir), width = 14, height = 16)
    print(p)
    dev.off()
    
    sel_tfs <- rownames(visium_slide)
    
    foreach(tf = sel_tfs) %dopar% {
        p <- SpatialFeaturePlot(object = visium_slide,
                                features = tf,
                                slot = "scale.data")
        pdf(sprintf("%s/hint_activity/HINT_%s.pdf", output_dir, tf))
        print(p)
        dev.off()
    }
    
    DefaultAssay(visium_slide) = 'dorothea'
    foreach(tf = sel_tfs) %dopar% {
        if(tf %in% rownames(visium_slide)){
            p <- SpatialFeaturePlot(object = visium_slide,
                                    features = tf)
            pdf(sprintf("%s/dorothea/Dorothea_Activity_%s.pdf", output_dir, tf))
            print(p)
            dev.off()
        }
    }
    
    DefaultAssay(visium_slide) = 'SCT'
    foreach(tf = sel_tfs) %dopar% {
        if(tf %in% rownames(visium_slide)){
            p <- SpatialFeaturePlot(object = visium_slide,
                                    features = tf)
            pdf(sprintf("%s/expression/SCT_%s.pdf", output_dir, tf))
            print(p)
            dev.off()
        }
    }

    if(!dir.exists(sprintf("%s/regulome", output_dir))){
        dir.create(sprintf("%s/regulome", output_dir))
    }    
    
    # check reguluom activity
    DefaultAssay(visium_slide) = 'SCT'
    
    # read TF-Target predicitons from Gene Target predictions of HINT
    dir_path = sprintf("../MotifToGene/%s/", sample)
    TF_networks = list()
    for (f in list.files(dir_path,"txt") ){
        name = sub(".txt","",f)
        name = sub(sample,"",name)
        df <- read.table(paste0(dir_path, f), header = FALSE)
        colnames(df) <- c("TF", "Gene", "Tag_Count")

        # only keep top 250 target genes for each  TF
        df <- df %>% group_by(TF) %>% top_n(n = 250, wt = Tag_Count)
        TF_networks[[name]] <- as.data.frame(df)
    }
    
    #transforming networks for dorothea input
    TF_networks_dor=list()
    for (n in names(TF_networks)){
        interactions = TF_networks[[n]]
        new_network = data.frame(tf=interactions[,1], 
                                 confidence = rep("A",dim(interactions)[1]),
                                 target=interactions[,2],
                                 mor=rep(1,dim(interactions)[1]))
        new_network = as.data.frame(sapply(new_network, toupper))
        new_network$mor = as.numeric(new_network$mor)  
        TF_networks_dor[[n]] = new_network 
    }
    
    ## We compute Viper Scores 
    Spatial_Regulomes_Cluster=list()
    for(n in names(TF_networks_dor)){
        interactions=TF_networks_dor[[n]]
        Spatial_Regulomes_Cluster[[n]] <- run_viper(as.matrix(visium_slide@assays$SCT@scale.data), 
                                                    interactions, 
                                                    options = list(method = "scale", minsize = 4, 
                                                                   eset.filter = FALSE, cores = 4, verbose = TRUE))
        visium_slide[[n]]= CreateAssayObject(counts = Spatial_Regulomes_Cluster[[n]])
        
        DefaultAssay(visium_slide) = n
        
        if(!dir.exists(sprintf("%s/regulome/%s", output_dir, n))){
            dir.create(sprintf("%s/regulome/%s", output_dir, n))
        }
        
        foreach(tf = rownames(visium_slide)) %dopar% {
            p <- SpatialFeaturePlot(object = visium_slide,
                                    features = tf)
            pdf(sprintf("%s/regulome/%s/%s.pdf", output_dir, n, tf))
            print(p)
            dev.off()
        }
        
    }
    saveRDS(visium_slide, sprintf("%s/%s.rds", output_dir, atac_to_spatial[sample]))
}
