# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' QC analysis of slides and scRNAseq: Is gene coverage similar?
#' 
#' 
#' Obtain expression profiles of slides and scell data sets
#' Pseudobulk them
#' Add them as tibbles
#' 


library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot)
library(stringr)
library(cowplot)
library(scater)
library(edgeR)
library(corrplot)


slides_ids = set_names(c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785"))

# CT-CT MISTy
visium_bulk = map(slides_ids, function(slide){
  
  print(slide)

  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file) 
  
  bulk_data = sumCountsAcrossCells(x = as.matrix(visium_slide@assays$SCT@counts),
                                   ids = visium_slide@meta.data[,"orig.ident"])
  
  bulk_df = tibble(id = slide, 
                   gene = rownames(bulk_data), 
                   expr_counts = bulk_data[,1])
  
  return(bulk_df)
})
 

scell_marker_files = list.files("./results/single_sample/all_markers")
scell_marker_files = unlist(lapply(strsplit(scell_marker_files,"\\."), 
                                   function(x) x[1]))

scell_dsets = set_names(scell_marker_files[-which(grepl("all",scell_marker_files))])

scell_summary = map(scell_dsets, function(scell_id){
  
  print(scell_id)
  
  scell_file = sprintf("./results/scell_RNA_final/%s.filtered.annotated.rds",
                       scell_id,scell_id)
  
  scell_data = readRDS(scell_file)
  
  bulk_data = sumCountsAcrossCells(x = as.matrix(scell_data@assays$RNA@counts),
                                   ids = scell_data@meta.data[,"orig.ident"])
  bulk_df = tibble(id = scell_id, 
                   gene = rownames(bulk_data), 
                   expr_counts = bulk_data[,1])
  
  variable_genes = scell_data@assays$RNA@var.features
  
  return(list("pseudobulk" = bulk_df,
              "var_genes" = variable_genes))
})

# Merging all data_sets

sample_dictionary = readRDS("./sample_dictionary.rds")

scell_bulk = lapply(scell_summary, function(x) x$pseudobulk)
  
all_expression = bind_rows(c(scell_bulk, visium_bulk)) %>%
  pivot_wider(id_cols = gene, 
              names_from = id,values_from = expr_counts)

all_expression_mat = as.matrix(all_expression[,-1])
rownames(all_expression_mat) = all_expression$gene

#Normalize as in Seurat
all_expression_mat[is.na(all_expression_mat)] <- 0
all_expression_mat = apply(all_expression_mat, 2, function(x) x/sum(x)) * 10000
all_expression_mat = log1p(all_expression_mat)

#Are there genes shared in all samples?
all_expression_mat_red = all_expression_mat[rowSums(all_expression_mat > 0)==15,]

#Scaling for PCA if needed
all_expression_scaled = scale(all_expression_mat_red)

cor_mat = cor(all_expression_mat_red, use = "pairwise")

corrplot(cor_mat,
         method = "color",
         tl.col = "black",
         is.corr = T)

#Reduced matrix
corrplot(cor_mat[grepl("1577",rownames(cor_mat)),1:7],
         method = "color",
         tl.col = "black",
         is.corr = F)

#Generate most variable gene filter
scell_vargenes = lapply(scell_summary, function(x) x$var_genes)

saveRDS(scell_vargenes %>% enframe(.), file = "./results/single_sample/all_vargenes.rds")











