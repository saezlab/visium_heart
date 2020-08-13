# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run misty in infarction,
#' I will use as features to predict,
#' differentially expressed genes from all clusters (top 10)
#' I will use ECM + cytokines for predictors

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

slides =  c("157779","157775")

for(slide in slides){
  print(slide)
  
  misty_ecm_out = sprintf("./results/single_slide/%s/misty/%s_ecm",
                          slide,slide)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_path_out = sprintf("./results/single_slide/%s/misty/%s_path",
                          slide,slide)
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  de_genes = dea_res[["SCT"]] %>%
    dplyr::filter(!grepl("MT-",gene)) %>%
    dplyr::arrange(cluster, p_val_adj, -avg_logFC) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice(1:20) %>%
    dplyr::select(gene) %>% pull() %>% unique()
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  # Run MISTy pipelines
  
  ls = c(2,5,10,20,50,100)
  
  # This is PROGENy high level analysis
  
  test_path = lapply(ls,para_ppln_seurat,
                     visium_slide = visium_slide,
                intra_assay = "SCT",
                intra_features = unique(de_genes),
                para_assay = "progeny",
                para_features = NULL,
                spot_ids = NULL,
                out_alias = misty_path_out)
  
  # This is ECM analysis
  
  # Prepare MISTy parameters
  hpredictors = readRDS("./markers/heart_predictors.rds")
  #Which ECM prots?
  slide_gex = visium_slide@assays$SCT@data
  slide_gex = as.matrix(slide_gex[hpredictors[hpredictors %in% rownames(slide_gex)],])
  predictor_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
  ligands = names(predictor_coverage[predictor_coverage>0.1])
  
  test_path = lapply(ls,para_ppln_seurat,visium_slide = visium_slide,
                    intra_assay = "SCT",
                    intra_features = unique(de_genes),
                    para_assay = "SCT",
                    para_features = ligands,
                    spot_ids = NULL,
                    out_alias = misty_ecm_out)
  
  # This is cytokine analysis
  
  # Prepare MISTy parameters
  gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")
  cytokines = gene_sets$MSIGDB_KEGG$KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
  hpredictors = cytokines
  slide_gex = visium_slide@assays$SCT@data
  slide_gex = as.matrix(slide_gex[hpredictors[hpredictors %in% rownames(slide_gex)],])
  predictor_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
  ligands = names(predictor_coverage[predictor_coverage>0.01])
  
  test_cytok = lapply(ls,para_ppln_seurat,visium_slide = visium_slide,
                    intra_assay = "SCT",
                    intra_features = unique(de_genes),
                    para_assay = "SCT",
                    para_features = ligands,
                    spot_ids = NULL,
                    out_alias = misty_cytok_out)
  
}

# Let's optimize everything

for(slide in slides){
  
  print(slide)
  
  ls = c(2,5,10,20,50,100)
  
  misty_ecm_out = sprintf("./results/single_slide/%s/misty/%s_ecm",
                          slide,slide)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_path_out = sprintf("./results/single_slide/%s/misty/%s_path",
                           slide,slide)
  
  print("path")
  get_optimal(out_dir_name = misty_path_out,
              ls = ls)
  
  print("ecm")
  get_optimal(out_dir_name = misty_ecm_out,
              ls = ls)
  
  print("cytok")
  get_optimal(out_dir_name = misty_cytok_out,
              ls = ls)
  
}

# Let's automatize summaries

for(slide in slides){
  
  print(slide)
  
  ls = c(2,5,10,20,50,100)
  
  misty_ecm_out = sprintf("./results/single_slide/%s/misty/%s_ecm",
                          slide,slide)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_path_out = sprintf("./results/single_slide/%s/misty/%s_path",
                           slide,slide)
  
  
  out_files = c(misty_ecm_out, misty_cytok_out, misty_path_out)
  
  for(f in out_files){
    
    MISTyout = MISTy_aggregator(results_folder = paste0(f,"_optim"),
                                             p.cutoff = 0.05)
    
    best_25 = MISTyout$performance %>% 
      dplyr::arrange(p.R2) %>% 
      dplyr::filter(p.R2 <= 0.15) %>% 
      dplyr::slice(1:25) %>%
      select(target) %>% pull()
    
    if(slide == "157775"){
    best_25 = unique(c(best_25, "CXCL8","C1QC","TNNI3",
                       "FN1","NCF1","ALOX5AP", "AREG",
                       "SAA1","PRKCB","BCL2A1","AKAP8"))
    }else{
    best_25 = unique(c(best_25, "EGFL7","SOX4","MYH6",
                       "TECRL","HERPUD1","ALOX5AP", "TAGLN",
                       "S100A9", "MT2A"))
    }
    
    plot_misty_performance(MISTy_out = MISTyout,
                           predicted_features = best_25,
                           pdf_out = paste0(f,
                                            "_optim_prfmnce.pdf"))
    
    plot_misty_importance(MISTy_out = MISTyout,
                          pdf_out = paste0(f,
                                           "_optim_imp.pdf"),
                          predicted_features = best_25,
                          predictors_features = NULL,
                          importance_cut = 1)
    
    
  }
  
}















