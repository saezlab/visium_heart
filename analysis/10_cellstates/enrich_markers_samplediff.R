# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we enrich differentially expressed genes 
#' between conditions to markers of states
#' 
#' The concept is quite simple, do genes associated to a 
#' disease group independent to the states?
#' 
#' 

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

annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic", "ischemic", "fibrotic"))

all_des <- all_des %>%
  left_join(annotation_names, by = c("groupA" = "patient_group")) %>%
  dplyr::select(-groupA) %>%
  rename("groupA" = patient_group_name) %>%
  left_join(annotation_names, by = c("groupB" = "patient_group")) %>%
  dplyr::select(-groupB) %>%
  rename("groupB" = patient_group_name)

# Upregulated genes -----------------

all_des_pos <- all_des %>% 
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

# Downregulated genes -----------------

all_des_neg <- all_des %>% 
  dplyr::filter(logFC < 0.5, FDR < 0.05) %>%
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

ora_analysis_pos <- all_des_pos %>%
  left_join(state_mrkrs) %>%
  dplyr::filter(!map_lgl(state_genes,is.null)) %>%
  dplyr::mutate(ora = map2(state_genes, cond_genes, function(state, cond){
    
    map(state, GSE_analysis, Annotation_DB = cond) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora) %>%
  dplyr::mutate(tested_genes = "upregulated")

ora_analysis_neg <- state_mrkrs %>%
  left_join(all_des_neg) %>%
  dplyr::mutate(ora = map2(state_genes, cond_genes, function(state, cond){
    
    map(state, GSE_analysis, Annotation_DB = cond) %>% 
      enframe() %>%
      unnest() %>%
      dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 
    
    
  })) %>%
  dplyr::select(cell_type, ora) %>%
  dplyr::mutate(tested_genes = "downregulated")


ora_plots_pos <-  ora_analysis_pos %>%
  mutate(dot_plts_up = map2(cell_type, ora, function(ct, o_dat) {
    
    o_dat <- o_dat %>%
      #dplyr::filter(corr_p_value < 0.1) %>%
      dplyr::mutate(log_p_value = -log10(corr_p_value))
      
    mapping_q <- quantile(o_dat$log_p_value, 0.95)
    
    o_dat <- o_dat %>%
      mutate(log_p_value = ifelse(log_p_value >= mapping_q, mapping_q, log_p_value))
    
    o_plt <- ggplot(o_dat, aes(y = gset, x = name, fill = log_p_value)) +
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
      coord_equal() +
      ggtitle("downregulated")
    
    
  })) %>%
  dplyr::select(dot_plts_up)

ora_plots_neg <-  ora_analysis_neg %>%
  mutate(dot_plts_neg = map2(cell_type, ora, function(ct, o_dat) {
    
    o_dat <- o_dat %>%
      #dplyr::filter(corr_p_value < 0.1) %>%
      dplyr::mutate(log_p_value = -log10(corr_p_value))
    
    mapping_q <- quantile(o_dat$log_p_value, 0.8)
    
    o_dat <- o_dat %>%
      mutate(log_p_value = ifelse(log_p_value >= mapping_q, mapping_q, log_p_value))
    
    
    o_plt <- ggplot(o_dat, aes(y = gset, x = name, fill = log_p_value)) +
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
      coord_equal() +
      ggtitle("upregulated")
    
    
  })) %>%
  dplyr::select(dot_plts_neg)

# Put them together and make a joint plot

all_plts <- left_join(ora_plots_pos, ora_plots_neg)

pdf("./results/cell_states/map_degs_all.pdf", height = 6, width = 5)

walk2(all_plts$dot_plts_up, all_plts$dot_plts_neg, function(a,b) {
  
  plot(plot_grid(a, b, align = "hv", ncol = 1))
  
  
})

dev.off()

rbind(ora_analysis_pos, ora_analysis_neg) %>%
  unnest() %>% write_csv("./results/cell_states/map_degs_all.csv")
