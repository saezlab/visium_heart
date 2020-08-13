# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

library(tidyverse)
library(Seurat)
library(clustree)
library(progeny)
library(dorothea)
library(cowplot)

source("visium_exploratory/MISTy_ct_runs.R")
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

# Just to get the spot IDs
get_CTcoordinates = function(ct_pattern, slide, max_filter = 0){
  ct_scores = ct_scores_lt_all[[slide]]
  ixs = grepl(ct_pattern,ct_scores$predicted_id,ignore.case = F)
  ct_scores = ct_scores[ixs,]
  ixs = ct_scores$prediction_score_max >= max_filter
  ct_ids = ct_scores$spot_id[ixs]
  return(ct_ids)
}

summarize_misty = function(slide,
                           ct_pattern,
                           optim_folder,
                           ct_scores,
                           p_value_trsh = 0.15,
                           R2_trsh = 0.1,
                           importance_cut = 1,
                           ident = "lt_id"){
  
  print(slide)
  print(ct_pattern)
  
  #Loading the slide
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  DefaultAssay(visium_slide) = "SCT"
  
  if(!is.null(ct_pattern)){
    
    #Get spots analyzed in MISTy for the selected ctpattern
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    #Reduce the object
    
    visium_slide = visium_slide[,spot_ids]
    
  }
  
  #Now let's get the features we want to analyze
  
  misty_out_folder = optim_folder
  
  #Let's define the pdf out
  
  group_pdf = paste0(misty_out_folder,"_group_imp.pdf")
  
  group_table = paste0(misty_out_folder,"_feature_ann.txt")
  
  #Then perform a differential expression analysis of the BIC genes
  
  Idents(visium_slide) = ident
  
  MISTyout = MISTy_aggregator(results_folder = misty_out_folder,
                              p.cutoff = 0.05)
  
  bic = MISTyout$performance %>% 
    dplyr::arrange(p.R2) %>% 
    dplyr::filter(p.R2 <= p_value_trsh,
                  intra.R2 >= R2_trsh) %>% 
    select(target) %>% pull()
  
  not_bic = MISTyout$targets[! MISTyout$targets %in% bic]
  
  if(length(bic)>0){
    
    annotation = FindAllMarkers(object = visium_slide,
                                test.use = "wilcox",
                                features = bic,
                                logfc.threshold = 0,
                                only.pos = T) %>%
      dplyr::select(gene,cluster)
    
    not_de = bic[!bic %in% annotation$gene]
    
    annotation = annotation %>%
      bind_rows(tibble(gene = not_de,
                       cluster = "shared"))
    
    colnames(annotation) = c("feature","ann")
    
    grouped_MISTy_out = group_MISTy_out(MISTy_out = MISTyout,
                    feature_ann = annotation)
    
    group_imp_summ = plot_misty_importance(MISTy_out = grouped_MISTy_out,
                          pdf_out = group_pdf,
                          importance_cut = importance_cut)
    
    write.table(annotation,
                file = group_table,
                quote = F,row.names = F,
                col.names = T,
                sep = "\t")
    
    return(group_imp_summ)
  } else{
    return(NULL)
  }
}

##################
# MAIN
###################

# 157781
# Cardiomyocytes
slide = "157781"
ct_pattern = "Cardiomyocytes"
p_dir = sprintf("./results/single_slide/%s/misty/",
                slide)
misty_dirs = list.dirs(p_dir)
misty_dirs = misty_dirs[grepl("optim",misty_dirs)]
misty_dirs = misty_dirs[grepl(ct_pattern,misty_dirs)]

for(dir in misty_dirs){
  print(dir)
  
  if(grepl("path",dir)){
    para_assay = "progeny"
  }else{
    para_assay = "SCT"
  }
  
  plot_misty_bic(MISTy_out_folders = dir,
                 p_value_trsh = 0.15,
                 R2_trsh = 0.05)
  
  g_imp = summarize_misty(slide = slide,ct_pattern = ct_pattern,
                  optim_folder = dir,ct_scores = ct_scores[[slide]],
                  p_value_trsh = 0.15,
                  R2_trsh = 0.05)
  
  # Plot para views
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  para_mat = get_para_matrix(visium_slide = visium_slide,
                             para_assay = para_assay,
                             para_features = unique(g_imp$para$Predictor),
                             l = 10)
  
  visium_slide[['paraview']] = CreateAssayObject(data = para_mat)
  DefaultAssay(visium_slide) = "paraview"
  
  pdf_file = paste0(dir,"_paraplts.pdf")
  
  pdf(file = pdf_file,width = 8,height = 7,onefile = T)
  
  for(P in unique(g_imp$para$Predictor)){
    
    plot(SpatialFeaturePlot(object = visium_slide,
                            features = P)) +
      theme(legend.position = "right")
    
  }
  
  dev.off()

}

# Fibroblasts
slide = "157781"
ct_pattern = "Fibroblasts"
p_dir = sprintf("./results/single_slide/%s/misty/",
                slide)
misty_dirs = list.dirs(p_dir)
misty_dirs = misty_dirs[grepl("optim",misty_dirs)]
misty_dirs = misty_dirs[grepl(ct_pattern,misty_dirs)]

for(dir in misty_dirs){
  print(dir)
  
  if(grepl("path",dir)){
    para_assay = "progeny"
  }else{
    para_assay = "SCT"
  }
  
  plot_misty_bic(MISTy_out_folders = dir,
                 p_value_trsh = 0.15,
                 R2_trsh = 0.05)
  
  g_imp = summarize_misty(slide = slide,ct_pattern = ct_pattern,
                  optim_folder = dir,ct_scores = ct_scores[[slide]],
                  p_value_trsh = 0.15,
                  R2_trsh = 0.05)
  
  # Plot para views
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  para_mat = get_para_matrix(visium_slide = visium_slide,
                             para_assay = para_assay,
                             para_features = unique(g_imp$para$Predictor),
                             l = 10)
  
  visium_slide[['paraview']] = CreateAssayObject(data = para_mat)
  DefaultAssay(visium_slide) = "paraview"
  
  pdf_file = paste0(dir,"_paraplts.pdf")
  
  pdf(file = pdf_file,width = 8,height = 7,onefile = T)
  
  for(P in unique(g_imp$para$Predictor)){
    
    plot(SpatialFeaturePlot(object = visium_slide,
                            features = P)) +
      theme(legend.position = "right")
    
  }
  
  dev.off()
}


# 157772

# Fibroblasts
slide = "157772"
ct_pattern = "Fibroblasts"
p_dir = sprintf("./results/single_slide/%s/misty/",
                slide)
misty_dirs = list.dirs(p_dir)
misty_dirs = misty_dirs[grepl("optim",misty_dirs)]
misty_dirs = misty_dirs[grepl(ct_pattern,misty_dirs)]

for(dir in misty_dirs){
  print(dir)
  
  if(grepl("path",dir)){
    para_assay = "progeny"
  }else{
    para_assay = "SCT"
  }
  
  
  plot_misty_bic(MISTy_out_folders = dir,
                 p_value_trsh = 0.15,
                 R2_trsh = 0.05)
  
  g_imp = summarize_misty(slide = slide,ct_pattern = ct_pattern,
                  optim_folder = dir,ct_scores = ct_scores[[slide]],
                  p_value_trsh = 0.15,
                  R2_trsh = 0.05)
  
  # Plot para views
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  para_mat = get_para_matrix(visium_slide = visium_slide,
                             para_assay = para_assay,
                             para_features = unique(g_imp$para$Predictor),
                             l = 10)
  
  visium_slide[['paraview']] = CreateAssayObject(data = para_mat)
  DefaultAssay(visium_slide) = "paraview"
  
  pdf_file = paste0(dir,"_paraplts.pdf")
  
  pdf(file = pdf_file,width = 8,height = 7,onefile = T)
  
  for(P in unique(g_imp$para$Predictor)){
    
    plot(SpatialFeaturePlot(object = visium_slide,
                            features = P)) +
      theme(legend.position = "right")
    
  }
  
  dev.off()
}


#"157779"
# Clusters
slide = "157779"
p_dir = sprintf("./results/single_slide/%s/misty/",
                slide)
misty_dirs = list.dirs(p_dir)
misty_dirs = misty_dirs[grepl("optim",misty_dirs)]

for(dir in misty_dirs){
  print(dir)
  
  if(grepl("path",dir)){
    para_assay = "progeny"
  }else{
    para_assay = "SCT"
  }
  
  plot_misty_bic(MISTy_out_folders = dir,
                 p_value_trsh = 0.15,
                 R2_trsh = 0.05)
  
  g_imp = summarize_misty(slide = slide,ct_pattern = NULL,
                  optim_folder = dir,
                  ct_scores = NULL,
                  p_value_trsh = 0.15,
                  R2_trsh = 0.05,
                  ident = "seurat_clusters")
  
  # Plot para views
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  para_mat = get_para_matrix(visium_slide = visium_slide,
                             para_assay = para_assay,
                             para_features = unique(g_imp$para$Predictor),
                             l = 10)
  
  visium_slide[['paraview']] = CreateAssayObject(data = para_mat)
  DefaultAssay(visium_slide) = "paraview"
  
  pdf_file = paste0(dir,"_paraplts.pdf")
  
  pdf(file = pdf_file,width = 8,height = 7,onefile = T)
  
  for(P in unique(g_imp$para$Predictor)){
    
    plot(SpatialFeaturePlot(object = visium_slide,
                            features = P)) +
      theme(legend.position = "right")
    
  }
  
  dev.off()
}

#"157775"
slide = "157775"

p_dir = sprintf("./results/single_slide/%s/misty/",
                slide)
misty_dirs = list.dirs(p_dir)
misty_dirs = misty_dirs[grepl("optim",misty_dirs)]

for(dir in misty_dirs){
  print(dir)
  
  if(grepl("path",dir)){
    para_assay = "progeny"
  }else{
    para_assay = "SCT"
  }
  
  plot_misty_bic(MISTy_out_folders = dir,
                 p_value_trsh = 0.15,
                 R2_trsh = 0.05)
  
  g_imp = summarize_misty(slide = slide,ct_pattern = NULL,
                  optim_folder = dir,ct_scores = ct_scores[[slide]],
                  p_value_trsh = 0.15,
                  R2_trsh = 0.05,
                  ident = "seurat_clusters")
  
  # Plot para views
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  para_mat = get_para_matrix(visium_slide = visium_slide,
                             para_assay = para_assay,
                             para_features = unique(g_imp$para$Predictor),
                             l = 10)
  
  visium_slide[['paraview']] = CreateAssayObject(data = para_mat)
  DefaultAssay(visium_slide) = "paraview"
  
  pdf_file = paste0(dir,"_paraplts.pdf")
  
  pdf(file = pdf_file,width = 8,height = 7,onefile = T)
  
  for(P in unique(g_imp$para$Predictor)){
    
    plot(SpatialFeaturePlot(object = visium_slide,
                            features = P)) +
      theme(legend.position = "right")
    
  }
  
  dev.off()
}
