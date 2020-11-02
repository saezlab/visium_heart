# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Add label transfer data to the seurat objects,
#' Add funcomics
#' Perform differential expression analysis
#' Assumes HINT activities were generated before

library(tidyverse)
library(Seurat)
source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

# Here I am trying to homogenize naming
ct_scores_lt_all = lapply(readRDS(file = "results/data_integration/integrated_meta.rds"),
                          function(x){
                            
                            x = x %>% rownames_to_column("spot_id") %>%
                              dplyr::mutate(predicted.id = gsub("[.]","_",
                                                                gsub(" ","_", predicted.id)))
                            colnames(x) = gsub("[.]","_", colnames(x))
                            
                            return(x)
                          })

slides_ids = set_names(c("157771", "157772", "157775",
                         "157777", "157779","157781",
                         "157782", "157785"))

for(slide in slides_ids){
  print(slide)
  
  slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                       slide,slide)
  
  dea_file = sprintf("./visium_results_manuscript/differential_analysis/%s_differential_feat.rds",
                     slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  #Adding funcomics - visiumtools implementation
  
  visium_slide = add_funcomics(visium_slide = visium_slide, 
                               top = 1000)
  
  #Adding cell_type scores from label transfer
  
  if(slide %in% slides_ids[names(ct_scores_lt_all)]){
    
    ct_scores = ct_scores_lt_all[[slide]]
    
    colnames(ct_scores) = gsub("prediction_score_","",
                               colnames(ct_scores)) 
    
    ct_scores_mat = as.matrix(ct_scores[,-c(1:5,ncol(ct_scores))])
    
    rownames(ct_scores_mat) = ct_scores[[1]]
    
    visium_slide[['labeltransfer']] = CreateAssayObject(counts = t(ct_scores_mat))
    
    #Performing differential expression analysis on known assays
    
    diff_feat = find_allfeat(visium_slide, assays_collection = c("SCT","dorothea","progeny",
                                                                 "labeltransfer","ECM"))
    
    # HINT requires lower LFC
    
    DefaultAssay(visium_slide) = "HINT_TF_Activity"
    
    daf_HINT = FindAllMarkers(visium_slide, 
                              logfc.threshold = 0.00005,
                              test.use = "wilcox")
    
    diff_feat$HINT_TF_Activity = daf_HINT
    
    DefaultAssay(visium_slide) = "SCT"
    
  } else{
    
    diff_feat = find_allfeat(visium_slide, assays_collection = c("SCT","dorothea",
                                                                 "progeny","ECM"))
    
  }
  
  saveRDS(visium_slide,file = slide_file)
  saveRDS(diff_feat, file = dea_file)
}











