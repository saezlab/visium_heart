# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Defining MISTy pipelines for all slides
#' 
#' 1) Tissue architecture: Are cells spatially constrained by others?
#' 2) CT ligands
#' 3) Pathway ligands
#' 
#' 

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(cowplot)
library(MISTy)
library(future)
source("./visium_exploratory/slide_processing.R")


run_ct_ct_MISTy = function(seurat_visium_obj,
                           l = 20,
                           out_dir_name = "default"){
  
  future::plan(multiprocess, workers = 6)
  
  # Getting data ready to create views
  
  geometry = seurat_visium_obj@images$slice1@coordinates
  
  # Data of cell-type scores, matched if possible
  if("ctscores_match" %in% Assays(seurat_visium_obj)){
    ct = as.matrix(seurat_visium_obj@assays$ctscores_match@data) %>% 
      t %>% data.frame 
  }else{
    ct = as.matrix(seurat_visium_obj@assays$ctscores@data) %>% 
      t %>% data.frame
  }
  
  # Order with geometry
  ct = ct[rownames(geometry),]
  
  # Change problematic labels
  colnames(ct) = gsub(pattern = "-",replacement = "_",
                      colnames(ct))
  
  # Creating views
  views_ct = create_initial_view(ct, unique.id = paste0(out_dir_name,"ct","_",l^2)) %>% 
    add_paraview(geometry[ ,2:3], l^2)
  
  MISTy_run = run_misty(views_ct,paste0(out_dir_name,"_",l^2))
  
}

plot_MISTy_ct = function(MISTy_out,out_pdf){
  
  # Supplemental 3. Improvement 
  R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
    mutate(target = targets) %>%
    pivot_longer(cols = -target, 
                 names_to = "name", 
                 values_to = "value")
  
  impr_plot = ggplot(R2_impr) +
    geom_point(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Change in variance explained") +
    xlab("Target") +
    #ylim(c(-5, 25)) +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  #Contribution
  
  coefs_plot = ggplot(MISTy_out$coefs) + 
    geom_col(aes(x=target, y=value, group=view, fill=view)) +
    xlab("Target") +
    ylab("Contribution") +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          legend.text = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     vjust = 0.5)) +
    scale_fill_brewer(palette="Dark2",
                      labels = c("Intrinsic","Para pathway" ,"Para ligands")) +
    labs(fill = "View")
  
  #Importance
  
  #Intra pathways
  importance_intra = tidyr::gather(MISTy_out$importance[[1]], "Predicted",
                                   "Importance", -Predictor)
  
  importance_intra_plt = ggplot(importance_intra,
                                aes(x = Predictor, 
                                    y = Predicted, 
                                    fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text = element_text(size=10),
          legend.key.size = unit(.6, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "white", 
                         mid = "white", 
                         high = scales::muted("blue"),
                         midpoint = 0.9) +
    xlab("Intrinsic pathways")
  
  # Examples from the paper
  
  #Para_pathways
  importance_para = tidyr::gather(MISTy_out$importance[[2]], "Predicted",
                                  "Importance", -Predictor)
  
  importance_para_plt = ggplot(importance_para,
                               aes(x = Predictor, 
                                   y = Predicted, 
                                   fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text = element_text(size=10),
          legend.key.size = unit(.6, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "white", 
                         mid = "white", 
                         high = scales::muted("blue"),
                         midpoint = 0.9) + 
    xlab("Para Pathways")
  
  pdf(file = out_pdf,width = 15,height = 15)
  plot(plot_grid(impr_plot,coefs_plot,
            importance_intra_plt, importance_para_plt,
            nrow = 2, ncol = 2))
  dev.off()
}


# Slides with current object
slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

# CT-CT MISTy
for(slide in slides_ids){
  
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  misty_out = sprintf("./results/single_slide/%s/mrun_ct_ct_%s",
                     slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  ls = c(2,5,10,20,50,100)
  
  test_A1 = lapply(ls, run_ct_ct_MISTy, seurat_visium_obj = visium_slide,
                   out_dir_name = misty_out)
  
  #get_optimal(out_dir_name = "./results/healthy_slides/mrun_ct_ct",ls = ls)
  
  #MISTy_out = MISTy_aggregator(results_folder = "./results/healthy_slides/mrun_ct_ct_optim")
   
}
  
# Get optimal, aggregate and plot

for(slide in slides_ids){
  
  print(slide)
  
  misty_out = sprintf("./results/single_slide/%s/mrun_ct_ct_%s",
                      slide,slide)
  
  misty_out_optim = sprintf("./results/single_slide/%s/mrun_ct_ct_%s_optim",
                      slide,slide)
  
  misty_out_pdf = sprintf("./results/single_slide/%s/mrun_ct_ct_%s.pdf",
                      slide,slide)
  
  ls = c(2,5,10,20,50,100)
  
  get_optimal(out_dir_name = misty_out,ls = ls)
  
  MISTy_out = MISTy_aggregator(results_folder = misty_out_optim)
  
  plot_MISTy_ct(MISTy_out = MISTy_out,out_pdf = misty_out_pdf)
  
}
























