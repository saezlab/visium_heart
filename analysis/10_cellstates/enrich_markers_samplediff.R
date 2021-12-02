# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we enrich differentially expressed genes 

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

# These are the markers of the differential expression analysis
  
all_des <- readRDS("./results/sample_comparison/all_cts_de_analysis.rds")

all_des <- all_des %>% 
  dplyr::filter(logFC > 0.5, FDR < 0.05) %>%
  dplyr::select(cell_type, groupA, groupB, gene) %>%
  mutate(id = paste0(groupA,"_vs_",groupB)) %>%
  dplyr::select(cell_type, gene, id) %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(data = map(data, function(dat){
    
    dat %>%
      group_by(id) %>%
      nest() %>%
      mutate(data = map(data, ~ .x[[1]])) %>%
      deframe()
    
  })) %>%
  dplyr::rename("cond_genes" = data)


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

ora_analysis <- state_mrkrs %>%
  left_join(all_des) %>%
  dplyr::mutate(ora = map2(state_genes, cond_genes, function(state, cond){
    
    map(state, GSE_analysis, Annotation_DB = cond) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora)


ora_plots <-  ora_analysis %>%
  mutate(dot_plts = map2(cell_type, ora, function(ct, o_dat) {
    
    o_dat <- o_dat %>%
      dplyr::mutate(log_p_value = -log10(corr_p_value)) %>%
      dplyr::mutate(log_p_value = ifelse(log_p_value == 0, NA, log_p_value))
      
    
    o_plt <- ggplot(o_dat, aes(y = gset, x = name, size = log_p_value)) +
      geom_point() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ylab("") +
      xlab("") +
      ggtitle(ct)
    
    
  })) %>%
  dplyr::select(dot_plts)


pdf("./results/cell_states/map_degs_all.pdf", height = 9, width = 9)

plot(cowplot::plot_grid(plotlist = ora_plots$dot_plts, align = "hv", nrow = 3))

dev.off()

ora_analysis %>% unnest() %>% write_csv("./results/cell_states/map_degs_all.csv")
