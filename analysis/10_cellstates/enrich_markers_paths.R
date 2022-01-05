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


gene_sets <- readRDS("./markers/Genesets_Dec19.rds")

hallmarks <- gene_sets$MSIGDB_HMARKS
canonical <- gene_sets$MSIGDB_CANONICAL


# Overrepresentation analysis

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

# Main

ora_analysis_hmarks <- state_mrkrs %>%
  dplyr::mutate(ora = map(state_genes, function(state){
    
    map(state, GSE_analysis, Annotation_DB = hallmarks) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora) %>%
  unnest() #%>%
  #dplyr::filter(corr_p_value < 0.2)

ora_analysis_canonical <- state_mrkrs %>%
  dplyr::mutate(ora = map(state_genes, function(state){
    
    map(state, GSE_analysis, Annotation_DB = canonical) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora) %>%
  unnest() #%>%
  #dplyr::filter(corr_p_value < 0.2)

# Per cell, arrange by name, GenesInList take top 10

useful_gsets_hmarks <- ora_analysis_hmarks %>%
  arrange(cell_type, name, -GenesInList) %>%
  group_by(name) %>%
  dplyr::slice(1:10) %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(gsets = map(data, function(dat) dat %>% pull(gset) %>% unique())) %>%
  dplyr::select(gsets) %>%
  dplyr::mutate(gsets = map(gsets, function(gst) {
    
    factor(gst, levels = gst)
    
  }))

useful_gsets_canonical <- ora_analysis_canonical %>%
  arrange(cell_type, name, -GenesInList) %>%
  group_by(name) %>%
  dplyr::slice(1:10) %>%
  group_by(cell_type) %>%
  nest() %>%
  dplyr::mutate(gsets = map(data, function(dat) dat %>% pull(gset) %>% unique())) %>%
  dplyr::select(gsets) %>%
  dplyr::mutate(gsets = map(gsets, function(gst) {
    
    factor(gst, levels = gst)
    
  }))

# Generate plots

ora_analysis_canonical_plts <- ora_analysis_canonical %>%
  nest() %>%
  left_join(useful_gsets_canonical) %>%
  dplyr::mutate(tile_plts = map2(data, gsets, function(dat, gst) {
    
    dat <- dat %>%
      dplyr::filter(gset %in% gst) %>%
      mutate(logpval = -log10(corr_p_value))
    
    mapping_q <- quantile(dat$logpval, 0.8)
    
    dat <- dat %>%
      mutate(logpval = ifelse(logpval >= mapping_q, mapping_q, logpval))
    
    dat %>%
      dplyr::filter(gset %in% gst) %>%
      ggplot(aes(y = factor(gset,
                            levels = levels(gst)),
                 x = name,
                 fill = logpval)) +
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
      ylab("") +
      xlab("") +
      coord_equal()
    
  }))

ora_analysis_hmarks_plts <- ora_analysis_hmarks %>%
  nest() %>%
  left_join(useful_gsets_hmarks) %>%
  dplyr::mutate(tile_plts = map2(data, gsets, function(dat, gst) {
    
    dat <- dat %>%
      dplyr::filter(gset %in% gst) %>%
      mutate(logpval = -log10(corr_p_value))
    
    mapping_q <- quantile(dat$logpval, 0.8)
    
    dat <- dat %>%
      mutate(logpval = ifelse(logpval >= mapping_q, mapping_q, logpval))
    
    dat %>%
      ggplot(aes(y = factor(gset,
                            levels = levels(gst)),
                 x = name,
                 fill = logpval)) +
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.text = element_text(size = 10), 
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
      ylab("") +
      xlab("") +
      coord_equal()
    
  }))

# Plots and files

ora_analysis_canonical %>%
  write_csv("./results/cell_states/ora_analysis_canonical_all.csv")

pdf("./results/cell_states/ora_analysis_canonical_all.pdf", height = 9, width = 11)

walk2(ora_analysis_canonical_plts$cell_type, ora_analysis_canonical_plts$tile_plts, function(ct, plt) {
  
  plot(plt + ggtitle(ct))
  
})

dev.off()


ora_analysis_hmarks %>%
  write_csv("./results/cell_states/ora_analysis_hmarks_all.csv")

pdf("./results/cell_states/ora_analysis_hmarks_all.pdf", height = 6, width = 11)

walk2(ora_analysis_hmarks_plts$cell_type, ora_analysis_hmarks_plts$tile_plts, function(ct, plt) {
  
  plot(plt + ggtitle(ct))
  
})

dev.off()








