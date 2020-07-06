# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Try to make zoom plots of slides
#' 

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)
library(cowplot)
library(scatterpie)
cols = c('Cardiomyocytes' = '#800000',
         'Cardiomyocytes 1' = '#800000',
         'Cardiomyocytes 2' = '#9A6324',
         'Cardiomyocytes 3' = '#808000',
         'Fibroblasts 1' = '#911eb4',
         'Fibroblasts 2' = '#e6beff',
         'Fibroblasts 3' = '#f032e6',
         'Fibroblasts 4' = '#f032e6',
         'Fibroblasts 5' = '#f032e6',
         'Fibroblasts 0' = '#f032e6',
         'Endothelial cells 1' = '#000075',
         'Endothelial cells 2' = 'blue',
         'Endothelial cells 3' = '#568198',
         'Endothelial cells 4' = '#469990',
         'Endothelial cells 5' = '#469990',
         'Endothelial cells 6' = '#469990',
         'Macrophages' = '#e6194B',
         'Macrophages 1' = '#e6194B',
         'Macrophages 2' = '#fabebe',
         'Pericytes' = '#f58231',
         'T cells' = '#ffe119',
         'Lymphatic endothelial cells' = '#ffd8b1',
         'Adipocytes' = '#000000',
         'Neuronal cells' = '#42d4f4',
         'Erythrocytes' = '#999999',
         'Endothelial cells (damaged)' = '#999999',
         'Vascular smooth muscle cells' = '#aaffc3')

names(cols) = gsub(" ","_", names(cols))
names(cols) = gsub("[(]","_", names(cols))
names(cols) = gsub("[)]","_", names(cols))


# Here I am trying to homogenize naming
ct_scores_lt_all = lapply(readRDS(file = "results/data_integration/integrated_meta.rds"),
                          function(x){
                            
                            x = x %>% rownames_to_column("spot_id") %>%
                              dplyr::mutate(predicted.id = gsub("[.]","_",
                                                                gsub(" ","_", predicted.id)))
                            colnames(x) = gsub("[.]","_", colnames(x))
                            
                            return(x)
                          })


# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")


plot_variable_spots = function(slide, 
                               ngenes = NULL,
                               vargenes = FALSE,
                               xlim,
                               ylim){
  require(scatterpie)
  require(Seurat)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  ct_scores = ct_scores_lt_all[[slide]]
  colnames(ct_scores) = gsub("prediction_score_","",colnames(ct_scores)) 

  ct_scores_mat = as.matrix(ct_scores[,-c(1:5,ncol(ct_scores))])
  
  rownames(ct_scores_mat) = ct_scores[[1]]
  
  if(vargenes){

  spot_ids = ct_scores[[1]]
    
  spot_var_ids = names(sort(apply(ct_scores_mat,1,var),
                  decreasing = F)[1:ngenes])
  
  # First generating new Ident
  id_data = set_names(colnames(visium_slide) %in% spot_var_ids,
                      colnames(visium_slide))
  
  visium_slide = AddMetaData(visium_slide,id_data,"id_data")
  Idents(visium_slide) = "id_data"
  
  }else{
    
    test_ids = visium_slide@images$slice1@coordinates[,c(2,3)] %>% dplyr::filter((row < xlim[2]) & (row>xlim[1]), 
                                                                                 (col > ylim[1]) & (col < ylim[2]))
    test_ids = rownames(test_ids)
    spot_ids = test_ids
    
    id_data = set_names(colnames(visium_slide) %in% spot_ids,
                        colnames(visium_slide))
    
    visium_slide = AddMetaData(visium_slide,id_data,"id_data")
    Idents(visium_slide) = "id_data"
    
    spot_var_ids = spot_ids
    
  }
  
  slide_plt = SpatialDimPlot(visium_slide, label = TRUE, 
                 label.size = 0,stroke = 0,label.box = F,alpha = 0) +
    scale_fill_manual(values = alpha(c("white","darkred"),0.4)) #+
    #theme(legend.position = "none")
  
  geometry = visium_slide@images$slice1@coordinates[,c(2,3)] %>%
    rownames_to_column("spot_id") %>%
    dplyr::filter(spot_id %in% spot_var_ids) %>%
    left_join(ct_scores[,-c(2,3,4,5,ncol(ct_scores))])
  
  set.seed(3997)
  
  ident_cols = sample(cols,ncol(ct_scores_mat))
  ident_cols = cols[colnames(geometry[,-c(1:3)])]
  
  pieplt = ggplot() + geom_scatterpie(aes(y=row, x=col, group=spot_id),pie_scale = 1.5, data=geometry,
                             cols=colnames(geometry)[-c(1,2,3)]) + coord_equal() +
    theme_minimal() + theme(legend.position = "bottom",
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                            panel.grid = element_blank()) + 
    scale_y_reverse() +
    scale_fill_manual(values = ident_cols)
  
  all_plt = cowplot::plot_grid(slide_plt,pieplt,ncol = 1,
                               rel_heights = c(0.5,1))
  
  return(all_plt)

}


xlim = c(20,35)
ylim = c(80,105)
slide = "157772"
slide_fibrosis = plot_variable_spots(slide = slide,
                                     ngenes = NULL,
                                     vargenes = FALSE,
                                     xlim = xlim,
                                     ylim = ylim)

pie_file = sprintf("./results/single_slide/%s/%s_pie.pdf",
                   slide,slide)

pdf(height = 11,width = 10, file = pie_file)

plot(slide_fibrosis)

dev.off()


xlim = c(25,40)
ylim = c(50,75)
slide = "157781"
slide_borderzone = plot_variable_spots(slide = slide,
                                     ngenes = NULL,
                                     vargenes = FALSE,
                                     xlim = xlim,
                                     ylim = ylim)
pie_file = sprintf("./results/single_slide/%s/%s_pie.pdf",
                   slide,slide)

pdf(height = 11,width = 10, file = pie_file)

plot(slide_borderzone)

dev.off()







