# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here, we will summarize the differential expression from differential_analysis.R
#' Performs overrepresentation analysis for DEA
#' We will generate supplemental material and figures 

library(tidyverse)
library(Seurat)
library(cowplot)
library(xlsx)
source("./visium/visiumtools/differential_test.R")
gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")

slides_ids = set_names(c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785"))

# Dot plots

for(slide in slides_ids){
  print(slide)
  
  dea_file = sprintf("./visium_results_manuscript/differential_analysis/%s_differential_feat.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  plot_file = sprintf("./visium_results_manuscript/dotplots/%s_de_features.pdf",
                      slide)
  
  DefaultAssay(visium_slide) = "SCT"
  
  # Dot plot of top 5 per cluster
  gene_summary = dea_res[["SCT"]] %>% 
    arrange(cluster,p_val_adj, avg_logFC) %>%
    dplyr::filter(avg_logFC > 0,
                  p_val_adj<=0.005) %>%
    dplyr::filter(!grepl("MT[-]",gene)) %>%
    dplyr::arrange(cluster,-avg_logFC) %>%
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
    arrange(cluster,-avg_logFC,p_val_adj) %>%
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
  possible_assays = c("progeny","ECM","labeltransfer")
  names(possible_assays) = possible_assays
  possible_assays = possible_assays[possible_assays %in%
                                      Assays(visium_slide)]
  
  plot_summary_data = lapply(dea_res[possible_assays], function(x){
    x %>% arrange(cluster,-avg_logFC,p_val_adj) %>%
      dplyr::filter(avg_logFC > 0,
                    p_val_adj<=0.005) %>%
      group_by(cluster) #%>%
      #slice(1:5)
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


# Now generate supplemental tables

for(slide in slides_ids){
  print(slide)
  
  slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                       slide,slide)
  
  dea_file = sprintf("./visium_results_manuscript/differential_analysis/%s_differential_feat.rds",
                     slide,slide)
  
  final_file = sprintf("./visium_results_manuscript/differential_analysis/%s_differential_features.xlsx",
                       slide,slide)
  
  dea_res = readRDS(dea_file)
  
  assays_res = names(dea_res)
  names(assays_res) = assays_res
  
  possible_assays = c("SCT","dorothea" ,"progeny",
                      "ECM","labeltransfer")
  
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
           p_val_adj < 0.0005) %>% 
    arrange(-avg_logFC,p_val_adj) %>%
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

  saveRDS(dea_res, file = dea_file)
}





