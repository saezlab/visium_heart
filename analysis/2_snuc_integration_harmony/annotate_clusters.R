# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we leverage the annotations of HCA and an external dataset to annotate lineages

library(Seurat)
library(tidyverse)
library(cowplot)

# Cluster markers from MI atlas
all_mrkrs <- readRDS("./processed_snrnaseq/integration/integrated_rnasamples_mrkrs.rds")[["RNA"]]

all_mrkrs <- all_mrkrs %>% 
  dplyr::select(cluster, avg_log2FC, p_val_adj, gene) %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  arrange(cluster, -avg_log2FC, p_val_adj, gene) %>%
  group_by(cluster) %>%
  dplyr::slice(1:300)

# Cluster markers from HCA atlas
hca_mrks <- read_csv("ext_data/hca_mrkrs.csv") %>%
  dplyr::filter(pvals_adj < 0.01) %>%
  arrange(group, -logfoldchanges) %>%
  dplyr::select(group, names) %>%
  group_by(group) %>%
  dplyr::slice(1:100) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]])) %>%
  deframe()

# Cluster markers from HF atlas

hf_kramann_atlas <- read_csv("ext_data/markers_HF.csv") %>%
  dplyr::rename("group" = cell_type) %>% dplyr::filter(pvals_adj < 0.01) %>%
  arrange(group, -logfoldchanges) %>%
  dplyr::select(group, names) %>%
  group_by(group) %>%
  dplyr::slice(1:100) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]])) %>%
  deframe()
  
# Now enrich with hypergeometric tests

# Function to do enrichment -----------------------------------------------------------------

GSE_analysis = function(geneList,Annotation_DB){
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  ResultsDF = ResultsDF %>% 
    rownames_to_column("gset") %>% 
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"), 
              as.numeric) %>% 
    dplyr::arrange(corr_p_value,GenesInList)
  
  return(ResultsDF)
  
}

# Annotate_clusters

mi_clust_anns_viahca <- all_mrkrs %>%
  dplyr::select(cluster, gene) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]])) %>%
  mutate(ora_res = map(data, GSE_analysis, Annotation_DB = hca_mrks)) %>%
  dplyr::select(ora_res) %>%
  unnest() %>%
  dplyr::select(cluster, gset, p_value)

hca_ann_plt <- mi_clust_anns_viahca %>% dplyr::slice(1:2) %>%
  ggplot(aes(x = cluster, y = gset, size = -log10(p_value))) +
  geom_point() +
  ylab("HCA_markers")

mi_clust_anns_viaHF <- all_mrkrs %>%
  dplyr::select(cluster, gene) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]])) %>%
  mutate(ora_res = map(data, GSE_analysis, Annotation_DB = hf_kramann_atlas)) %>%
  dplyr::select(ora_res) %>%
  unnest() %>%
  dplyr::select(cluster, gset, p_value)

hf_ann_plt <- mi_clust_anns_viaHF %>% dplyr::slice(1:2) %>%
  ggplot(aes(x = cluster, y = gset, size = -log10(p_value))) +
  geom_point() +
  ylab("HF_atlas_markers")

annotations_plt <- plot_grid(hca_ann_plt, hf_ann_plt, 
                         nrow = 2, align = "hv")

pdf(file = "processed_snrnaseq/integration/atlas_annotation.pdf",
    width = 10,
    height = 4.5)

plot(annotations_plt)

dev.off()

write.table(all_mrkrs, quote = F,sep = "\t",file = "processed_snrnaseq/integration/clstr_mrkrs.txt",row.names = F, col.names = T)

# Facilitate annotations:

mi_proposal_viahf <- mi_clust_anns_viaHF %>% 
  dplyr::slice(1:2) %>% 
  dplyr::filter(p_value < .0001) %>%
  dplyr::rename("opt_clust_integrated" = cluster,
         "cell_type" = gset) %>%
  dplyr::select(-p_value) %>%
  write.table(., 
              quote = F, 
              sep = "\t",
              file = "processed_snrnaseq/integration/annotations_proposal_viahf.txt",
              row.names = F, 
              col.names = T)
  
mi_proposal_viahca <- mi_clust_anns_viahca %>% 
  dplyr::slice(1:2) %>% 
  dplyr::filter(p_value < .0001) %>%
  dplyr::rename("opt_clust_integrated" = cluster,
                "cell_type" = gset) %>%
  dplyr::select(-p_value) %>%
  write.table(., 
              quote = F, 
              sep = "\t",
              file = "processed_snrnaseq/integration/annotations_proposal_viahca.txt",
              row.names = F, 
              col.names = T)
