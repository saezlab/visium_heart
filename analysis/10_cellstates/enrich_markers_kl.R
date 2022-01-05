# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we enrich differentially expressed genes 
#' between conditions to a collection of canonical pathways

library(tidyverse)
library(cowplot)

# These are the markers of the cell-states

state_mrkrs <- tibble(marker_file = list.files("./cell_states", full.names = T)) %>%
  dplyr::mutate(cell_type = gsub("./cell_states/", "", marker_file)) %>%
  dplyr::mutate(marker_file = paste0(marker_file,"/annotation.rds")) %>%
  dplyr::mutate(markers = map(marker_file, readRDS)) %>%
  dplyr::select(cell_type, markers) %>%
  unnest() %>%
  dplyr::select(cell_type, p_val_adj, cluster, gene) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(-p_val_adj) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(data = map(data, function(dat){
    
    dat %>%
      group_by(cluster) %>%
      nest() %>%
      mutate(data = map(data, ~ .x[[1]])) %>%
      deframe()
    
    
  })) %>%
  dplyr::rename("state_genes" = data)

# Then we need the markers of KL atlas

kl_annotation <- read_csv("markers/kl_miatlas_clustermrkrs_KL.csv") %>%
  na.omit() %>%
  dplyr::select(cluster, annotation) %>%
  dplyr::mutate(annotation = gsub(" ","_", annotation) %>%
                  gsub(":","_", .) %>%
                  gsub("-","_", .))

kl_markers <- read_csv("markers/kl_miatlas_clustermrkrs_KL.csv") %>%
  dplyr::select(-annotation) %>%
  left_join(kl_annotation) %>%
  dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  dplyr::select(gene, annotation) %>%
  group_by(annotation) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

# Enrichment analysis

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

###

ora_analysis_kl <- state_mrkrs %>%
  dplyr::mutate(ora = map(state_genes, function(state){
    
    map(state, GSE_analysis, Annotation_DB = kl_markers) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora) %>%
  unnest() %>%
  dplyr::filter(corr_p_value < 0.05) %>%
  dplyr::mutate(label = paste0(cell_type, "_", name)) %>% 
  group_by(cell_type) %>%
  nest()

## Get the top 2 annotations per label

ora_analysis_kl_filt <- ora_analysis_kl %>%
  dplyr::mutate(data = map(data, function(dat) {
    
    
    dat %>%
      arrange(label, corr_p_value) %>%
      group_by(label) %>%
      slice(1:2)
    
    
    
  }))

enrichment_plots <- map2(ora_analysis_kl_filt$cell_type, ora_analysis_kl_filt$data, function(ct, dat) {
  
  dat %>%
    ggplot(aes(x = label, y = gset, size = -log10(corr_p_value))) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5)) +
    ggtitle(ct) +
    xlab("MI cell-types") +
    ylab("KL cell-types") +
    coord_equal()
  
})

pdf("./results/cell_states/ora_analysis_KL.pdf", height = 4, width = 5)

walk(enrichment_plots, plot)

dev.off()

#

ora_analysis_kl %>% unnest() %>% write_csv(., "./results/cell_states/ora_analysis_KL.csv")



