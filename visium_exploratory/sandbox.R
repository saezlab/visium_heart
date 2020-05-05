library(Seurat)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(purrr)
library(cowplot)
library(clustree)
library(xlsx)
library(cowplot)
source("./visium_exploratory/slide_processing.R")

ct_scores_lt_all = readRDS(file = "results/data_integration/rna_meta.rds")

# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

# Quality plot
slide = "157775"

for(slide in slides_ids){
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  SpatialFeaturePlot(visium_slide, feature = c("nFeature_Spatial",
                                               "nCount_Spatial"),ncol = 1)
  
  SpatialFeaturePlot(visium_slide, feature = c("nFeature_Spatial"),ncol = 1) + 
    theme(legend.position = "right")
  
  
  
  
}


ct_scores_lt_all = readRDS(file = "results/data_integration/rna_meta.rds")


slide = "157771"
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)

visium_meta = visium_slide@meta.data %>% rownames_to_column("spot_id")

ct_scores_lt = ct_scores_lt_all[[slide]]
ct_scores_ftprnt = t(as.matrix(visium_slide@assays$ctscores@data))

ct_scores_ftprnt = pnorm(ct_scores_ftprnt[rownames(ct_scores_lt),])

ct_scores_ftprnt_match = t(as.matrix(visium_slide@assays$ctscores_match@data))
ct_scores_ftprnt_match = pnorm(ct_scores_ftprnt_match[rownames(ct_scores_lt),])

# Exploration

plot(ct_scores_lt$prediction.score.max,
     ct_scores_lt$nCount_SCT,cex=0.3)

plot(apply(ct_scores_ftprnt_match,1,max)[visium_meta$spot_id], 
     visium_meta$nCount_SCT,cex=0.3)

plot(apply(ct_scores_ftprnt_match,1,max)[rownames(ct_scores_lt)], 
     ct_scores_lt$prediction.score.max,cex=0.3)



markers_mat = readRDS("./markers/markers_matrix.rds")



## Association of Pathways,ct-scores and TFs

slide = "157772"
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)

ct_data = t(as.matrix(visium_slide@assays$ctscores@data)) %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>% arrange(spot_id)

tf_data = t(as.matrix(visium_slide@assays$dorothea@data)) %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>% arrange(spot_id)

path_data = t(as.matrix(visium_slide@assays$progeny@data)) %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>% arrange(spot_id)

all_tfs = colnames(tf_data)[-1]
names(all_tfs) = all_tfs

TGFB_res = data.frame(t(sapply(all_tfs, function(tf){
  
  model_df = data_frame(path = path_data$TGFb,
                        TF = tf_data[,tf])
  
  model_res = lm(path~TF,data = model_df)
  lmres = summary(model_res)$coefficients[2,c(3,4)]
  names(lmres) = c("t","p_value")

  return(lmres)
}
))) %>% rownames_to_column("TF") %>%
  arrange(p_value)

TGFB_associatedTFs = c("SMAD1","SMAD4","ID1",
                       "SMAD6","TGIF2","E2F4",
                       "TFDP1","TGIF2","PITX2")

TGFB_res %>% dplyr::filter(TF %in% TGFB_associatedTFs)

# Fibroblast_associated: ARNT

#Cell type scores
DefaultAssay(visium_slide) = "ctscores"
fibro_plot = SpatialFeaturePlot(object = visium_slide,
                        features = "fibroblasts")

#PROGENy scores
DefaultAssay(visium_slide) = "progeny"
tgfb_plot = SpatialFeaturePlot(object = visium_slide,
                        features = "TGFb")

#TF scores
DefaultAssay(visium_slide) = "dorothea"
TF_plot = SpatialFeaturePlot(object = visium_slide,
                               features = c("E2F4","ARNT"))

plot_grid(plot_grid(fibro_plot,tgfb_plot,nrow = 1),
          TF_plot,ncol = 1)


# TF regulons
dorothea_regulon_human = read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")
# We obtain the regulons based on interactions with confidence level A, B and C
regulon = dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C")) %>%
  df2regulon()

TGFb_regulons = regulon[names(regulon) %in% TGFB_associatedTFs]
TGFb_regulons =lapply(TGFb_regulons, function(x) names(x$tfmode))
E2F4_regulon = TGFb_regulons$E2F4

# Know back to genes
moran_res = read.csv("./results/spatial_cor/moran_intersections.csv", 
                     header = T,
                     stringsAsFactors = F) %>% 
  dplyr::filter(Group == "Chronic")

spark_res = read.csv("./results/spatial_cor/spark_intersections.csv", 
                     header = T,
                     stringsAsFactors = F) %>% dplyr::filter(Group == "Chronic") 

spark_res %>% dplyr::filter(Gene %in% unlist(TGFb_regulons))
spark_res %>% dplyr::filter(Gene %in% E2F4_regulon &
                            ! Gene %in% moran_res$Gene) %>% 
  dplyr::select(Gene) %>% pull()

DefaultAssay(visium_slide) = "SCT"
gene_plot = SpatialFeaturePlot(object = visium_slide,
                             features = c("JUN","PPIA"))

gene_plot = SpatialFeaturePlot(object = visium_slide,
                               features = c("HMGN2","PPIA"))

gene_plot = SpatialFeaturePlot(object = visium_slide,
                               features = c("HMGB1","PPIA"))



plot_grid(plot_grid(fibro_plot,tgfb_plot,nrow = 1),
          TF_plot,ncol = 1)


spark_res$Gene

SpatialFeaturePlot(object = visium_slide,
                    features =c("nFeature_Spatial","nCount_Spatial"),)



#Slide: 157785


slide = "157785"
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide_validation = readRDS(slide_file)

#Cell type scores
DefaultAssay(visium_slide_validation) = "ctscores"
fibro_plot = SpatialFeaturePlot(object = visium_slide_validation,
                                features = "fibroblasts")

#PROGENy scores
DefaultAssay(visium_slide_validation) = "progeny"
tgfb_plot = SpatialFeaturePlot(object = visium_slide_validation,
                               features = "TGFb")

#TF scores
DefaultAssay(visium_slide_validation) = "dorothea"
TF_plot = SpatialFeaturePlot(object = visium_slide_validation,
                             features = c("E2F4","ARNT"))

upper_val = plot_grid(plot_grid(fibro_plot,tgfb_plot,nrow = 1),
          TF_plot,ncol = 1)

######
DefaultAssay(visium_slide_validation) = "SCT"


gene_plot_val = SpatialFeaturePlot(object = visium_slide_validation,
                               features = c("HMGN2","PPIA"))

gene_plot_val = SpatialFeaturePlot(object = visium_slide_validation,
                               features = c("HMGB1","PPIA"))

plot_grid(upper_val,
          gene_plot_val,ncol = 1)

######


sample_dictionary = readRDS("./sample_dictionary.rds")
ct_scores_lt_all = readRDS(file = "results/data_integration/integrated_meta.rds")

# Healthy example

slide = "157771"
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)

ct_scores_lt = ct_scores_lt_all[[slide]][,7:15]
colnames(ct_scores_lt) = gsub("prediction.score.","",
                              colnames(ct_scores_lt))

ct_scores_lt = t(as.matrix(ct_scores_lt[colnames(visium_slide),]))

visium_slide[['ct_scores_lt']] = CreateAssayObject(counts = ct_scores_lt)


DefaultAssay(visium_slide) = "ct_scores_lt"

SpatialFeaturePlot(visium_slide, features = "Vascular.smooth.muscle.cells")


# Other Fibrotic

slide = "157785"
slide_file = sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)

ct_scores_lt = ct_scores_lt_all[[slide]][,7:13]
colnames(ct_scores_lt) = gsub("prediction.score.","",
                              colnames(ct_scores_lt))

ct_scores_lt = t(as.matrix(ct_scores_lt[colnames(visium_slide),]))

visium_slide[['ct_scores_lt']] = CreateAssayObject(counts = ct_scores_lt)


DefaultAssay(visium_slide) = "ct_scores_lt"

SpatialFeaturePlot(visium_slide, features = "fibroblast")


# Quality DF =

slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

names(slides_ids) = slides_ids

QC_basic_stats = lapply(slides_ids, function(slide){
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  return(visium_slide@meta.data[,c("nCount_Spatial","nFeature_Spatial")])
  
})

AllMeta_data = bind_rows(QC_basic_stats,.id = "slide") %>%
  group_by(slide) %>% mutate("n_spots" = length(nCount_Spatial))

df = bind_rows(QC_basic_stats,.id = "slide")
df.plot1 = df %>% group_by(slide) %>% summarise(num_spot = n())

p1 <- ggplot(data = df.plot1, aes(x = slide, y = log10(num_spot))) +
  geom_bar(aes(color = slide, fill = slide), 
           alpha = 1, stat = "identity") +
  xlab("") + ylab("Number of \nspots (log10)") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p2 <- ggplot(data = df, aes(x = slide, y = log10(nCount_Spatial))) +
  geom_boxplot(aes(color = slide)) +
  geom_violin(aes(color = slide, fill = slide), alpha = 0.5) +
  xlab("") + ylab("Counts (log10)") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p3 <- ggplot(data = df, aes(x = slide, y =  log10(nFeature_Spatial))) +
  geom_boxplot(aes(color = slide)) +
  geom_violin(aes(color = slide, fill = slide), alpha = 0.5) +
  xlab("") + ylab("Genes (log10)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none")


p <- plot_grid(p1, p2, p3, ncol = 1, align = "v", 
               rel_heights = c(1, 1, 1.3))

pdf("./results/qc_Spatial.pdf", width = 4, height = 6)
print(p)
dev.off()

### Fix ORA mistake

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

for(slide in slides_ids){
  print(slide)
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  final_file = sprintf("./results/single_slide/%s/%s_ora_fixed.xlsx",
                       slide,slide)

  dea_res = readRDS(dea_file)
  
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










