# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Defining MISTy pipelines for all slides for specific cell-types
#' 
#' 1) First, identify cell types in space based on label transfer
#' 2) From those IDs, fit an intra view with marker genes that overlap (n=25)
#' 3) Fit a paraview with ligands or pathways (parameter) and extract the ids
#' from the intraview
#' 4) Run MISTy
#' 
#' For each slide, we need assignment and marker genes 
#' 
#' 
#' 
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
library(MISTy)
library(future)
source("./visium_exploratory/slide_processing.R")

scell_marker_files = list.files("./results/single_sample/all_markers")
scell_marker_files = unlist(lapply(strsplit(scell_marker_files,"\\."), 
                                   function(x) x[1]))
sample_dictionary = readRDS("./sample_dictionary.rds")
all_var_genes = readRDS(file = "./results/single_sample/all_vargenes.rds") %>% 
  enframe(.,"id") %>% unnest()

# Here I am trying to homogenize naming
ct_scores_lt_all = lapply(readRDS(file = "results/data_integration/integrated_meta.rds"),
                          function(x){
                            
                            x = x %>% rownames_to_column("spot_id") %>%
                              dplyr::mutate(predicted.id = gsub("[.]","_",
                                                                gsub(" ","_", predicted.id)))
                            colnames(x) = gsub("[.]","_", colnames(x))
                            
                            return(x)
                          })

## Function to get MISTy results
## Inputs:
## seurat_visium_obj: visium_slide
## spot_ids: selected spots considered in MISTy
## intrinsic_markers: cell-type markers that define the intra view
## ligands: if NULL, para view built with pathways, vectro
## pathways: if NULL, para view built with ligands,
## l = radius parameter
## out_dir_name = for misty out
## Outputs:
## MISTy folders
run_mrkr_MISTy = function(seurat_visium_obj,
                            spot_ids,
                            intrinsic_markers,
                            ligands = NULL,
                            pathways = NULL,
                            l = 20,
                            out_dir_name = "default"){
  
  plan(multiprocess, workers = 6)
  
  # Getting data ready to create views
  geometry = seurat_visium_obj@images$slice1@coordinates[,c(2,3)]
  
  # Intrinsic main view: Marker genes of a slide
  markers_gex = as.matrix(seurat_visium_obj@assays$SCT@data)
  intrinsic_markers = intrinsic_markers[intrinsic_markers %in% rownames(markers_gex)]
  
  markers_gex = markers_gex[intrinsic_markers,spot_ids] %>%
    t %>% data.frame(check.names = F)
  
  views_main = create_initial_view(markers_gex, unique.id = "intra")
  
  if(length(pathways)>0){
  
    pth = as.matrix(seurat_visium_obj@assays$progeny@data) %>% 
      t %>% data.frame(check.names = F)
    
    pth = pth[rownames(geometry),pathways]
    
    views_para = create_initial_view(pth, unique.id = paste0(out_dir_name,"path","_",l^2)) %>% 
      add_paraview(geometry[ ,1:2], l^2)
    
    # Final view = comes from the view above
    data_red = views_para[[3]]$data
    rownames(data_red) = rownames(pth) #we named rows just for easy access
    data_red = data_red[rownames(markers_gex),]
    
  }else if(length(ligands)>0){
    
    ligs = as.matrix(seurat_visium_obj@assays$SCT@data) %>% 
      t %>% data.frame(check.names = F)
    
    ligands = ligands[ligands %in% colnames(ligs)]
    
    ligs = ligs[rownames(geometry),ligands]
    
    views_para = create_initial_view(ligs, unique.id = paste0(out_dir_name,"ligs","_",l^2)) %>% 
      add_paraview(geometry[ ,1:2], l^2)
    
    # Final view = comes from the view above
    data_red = views_para[[3]]$data
    rownames(data_red) = rownames(ligs) #we named rows just for easy access
    data_red = data_red[rownames(markers_gex),]
    
  }
  
  views = views_main %>% 
    add_views(create_view(paste0("para_",l^2),
                          data_red))
  
  MISTy_run = run_misty(views,paste0(out_dir_name,"_",l^2))
  
}

#
# Identifies IDs of spots that belong to the ct_pattern
get_CTcoordinates = function(ct_pattern, slide, max_filter = 0){
  ct_scores = ct_scores_lt_all[[slide]]
  ixs = grepl(ct_pattern,ct_scores$predicted_id,ignore.case = F)
  ct_scores = ct_scores[ixs,]
  ixs = ct_scores$prediction_score_max >= max_filter
  ct_ids = ct_scores$spot_id[ixs]
  return(ct_ids)
}

#Plots MISTy output
plot_MISTy_gene_path = function(all_markers_ann,
                                MISTy_out,
                                pdf_out,
                                ligs,
                                performance){
  #Generating figures
  
  # First colors and order
  set.seed(30997)
  color_df = unique(all_markers_ann$cluster)
  color_v = sample((RColorBrewer::brewer.pal(n = 8,
                                     "Dark2")), length(color_df))
  
  color_df = tibble(cluster = color_df, ccolor = color_v)
  color_df = left_join(all_markers_ann, color_df)
  
  color_df = color_df %>% mutate(ccolor = ifelse(shared==T,"black",
                                                 ccolor),
                                 cluster = ifelse(shared==T,"shared",
                                                  cluster)) %>%
    dplyr::select(gene,cluster,ccolor) %>% unique()
  
  # Overall performance
  
  performance_plt = performance %>% tidyr::pivot_longer(cols = -target) %>%
    dplyr::filter(grepl("R2",name)) %>%
    dplyr::mutate(target = factor(target,
                                  levels = color_df$gene)) %>%
    ggplot(aes(fill = name, y = target, x = value)) +
    geom_bar(stat = "identity",position="dodge") + 
    theme_minimal() +
    theme(axis.text.y = element_text(color = color_df$ccolor))
  
  R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
    mutate(target = targets) %>%
    pivot_longer(cols = -target, 
                 names_to = "name", 
                 values_to = "value") %>% arrange(desc(value)) %>%
    mutate(target = factor(target,
                           levels = color_df$gene))
  
  impr_plot = ggplot(R2_impr) +
    geom_point(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Change in variance explained") +
    xlab("Target") +
    #ylim(c(-5, 25)) +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1,
                                     color = color_df$ccolor))
  
  #Look at importances of views
  contribution_df = MISTy_out$coefs %>% dplyr::filter(target %in% R2_impr$target)
  coefs_plot = ggplot(MISTy_out$coefs) + 
    geom_col(aes(x=factor(target, levels = color_df$gene), 
                 y=value, group=view, fill=view)) +
    xlab("Target") +
    ylab("Contribution") +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          legend.text = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     vjust = 0.5)) +
    scale_fill_brewer(palette="Dark2",
                      labels = c("Intrinsic","Para pathway")) +
    labs(fill = "View")
  
  
  #Intra
  importance_intra = tidyr::gather(MISTy_out$importance[[1]], "Predicted",
                                   "Importance", -Predictor) %>%
    left_join(all_markers_ann, by = c("Predicted"="gene")) %>%
    mutate(PredictedCell = cluster) %>% 
    select(Predictor,Predicted,
           Importance,PredictedCell) %>%
    left_join(all_markers_ann, by = c("Predictor"="gene")) %>%
    mutate(PredictorCell = cluster) %>% 
    select(Predictor,Predicted,Importance,
           PredictorCell,
           PredictedCell) %>%
    arrange(PredictedCell,PredictorCell,Predicted,Predictor)
  
  importance_intra = importance_intra %>%
    mutate(Predictor = factor(Predictor,
                              levels = color_df$gene),
           Predicted = factor(Predicted,
                              levels = color_df$gene))
  
  importance_intra_plt = ggplot(importance_intra,
                                aes(x = Predictor, 
                                    y = Predicted, 
                                    fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     color = color_df$ccolor),
          axis.text.y = element_text(color = color_df$ccolor),
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
  
  #Para
  
  importance_para = tidyr::gather(MISTy_out$importance[[2]], "Predicted",
                                  "Importance", -Predictor) %>%
    dplyr::mutate(Predicted = factor(Predicted,
                                     levels = color_df$gene)) %>%
    arrange(-Importance)
  
  if(ligs == TRUE){
    
    useful_predictors = importance_para %>% 
      arrange(Predicted,-Importance) %>% 
      dplyr::group_by(Predicted) %>% 
      dplyr::slice(1:10) %>% 
      dplyr::select("Predictor") %>%
      dplyr::pull() %>%
      unique()
    
    importance_para = importance_para %>%
      dplyr::filter(Predictor %in% useful_predictors) %>%
      arrange(-Importance)
    
    write.table(importance_para, col.names = T,
                row.names = F,sep = ",",
                file = gsub("pdf","csv",pdf_out))
    
    ligs_mat = MISTy_out$importance[[2]] %>% 
      dplyr::filter(Predictor %in% useful_predictors)
    
    ligs_mat_names = ligs_mat$Predictor

    ligs_mat = as.matrix(ligs_mat[,-which(grepl("Predictor",colnames(ligs_mat)))])
    
    rownames(ligs_mat) = ligs_mat_names
    
    ligs_clust =  hclust(d = t(dist(ligs_mat)))
    
    ligands_order = ligs_clust$labels[ligs_clust$order]
  
    importance_para_plt = ggplot(importance_para,
                                 aes(x = factor(Predictor,
                                                levels = ligands_order), 
                                     y = Predicted, 
                                     fill = Importance)) + geom_tile() + 
      theme(panel.grid = element_blank(),
            axis.title = element_text(size=11),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.text.y = element_text(color = color_df$ccolor),
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
  }else{
    
  write.table(importance_para, col.names = T,
                row.names = F,sep = ",",
                file = gsub("pdf","csv",pdf_out))
  
  importance_para_plt = ggplot(importance_para,
                               aes(x = Predictor, 
                                   y = Predicted, 
                                   fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text.y = element_text(color = color_df$ccolor),
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
  }
  
  pdf(file = pdf_out, width = 12, height = 9,onefile = TRUE)
  plot(performance_plt)
  plot(impr_plot)
  plot(coefs_plot)
  plot(importance_intra_plt)
  plot(importance_para_plt)
  
  dev.off()
  
}

###########################
#
# MAIN FUNCTION
#
###########################

# CT-CT MISTy
MISTy_marker_pipeline = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  
  #Defining out files
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  misty_out = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_path_%s",
                      slide,slide,ct_pattern)
  
  misty_out_ligs = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_ligs_%s",
                      slide,slide,ct_pattern)
  
  QC_out = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_path_%s_qc.pdf",
                    slide,slide,ct_pattern)
  
  pdf_out = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_path_%s.pdf",
                    slide,slide,ct_pattern)
  
  pdf_out_ligs = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_ligs_%s.pdf",
                         slide,slide,ct_pattern)
  
  # General parameters
  
  visium_slide = readRDS(slide_file)
  
  ls = c(2,5,10,20,50,100)
  
  #Cell type scores from specific scRNA if data available
  scell_id = dplyr::filter(sample_dictionary, visium_ids == slide) %>%
    dplyr::select(scell_ids) %>% pull()
  
  if(scell_id %in% scell_marker_files){
    
    #Reading genesorter markers
    scell_markers = read.table(file = sprintf("./results/single_sample/top50/%s.top50.genes.txt",
                                              scell_id),
                               header = T, sep = "\t",
                               stringsAsFactors = F) %>%
      tidyr::pivot_longer(cols=colnames(.),"cluster",values_to = "gene") %>%
      dplyr::mutate(cluster = gsub("[.]","_",cluster))
    
    
    # # Here create score matrix
    #  scell_markers = read.table(file = sprintf("./results/single_sample/all_markers/%s.AllMarkers.txt",
    #                                            scell_id),
    #                             header = T, sep = "\t",
    #                             stringsAsFactors = F) %>%
    #    dplyr::mutate(cluster = gsub(" ","_",cluster)) %>%
    #    dplyr::mutate(cluster = gsub("\\+","pos",cluster)) %>%
    #    group_by(cluster) %>% arrange(p_val_adj) %>%
    #    dplyr::slice(1:200) %>% dplyr::filter(!grepl("-",gene))
    
    # Filter by most variable genes in the specific slide
    spec_vargns = all_var_genes %>% dplyr::filter(id == scell_id) %>%
      dplyr::select(value) %>% pull()
    
    scell_markers = dplyr::filter(scell_markers, gene %in% spec_vargns)
    
    # Continue processing of markers
    available_cts = unique(scell_markers$cluster)
    
    useful_cts = available_cts[grepl(pattern = ct_pattern, available_cts,
                       ignore.case = F)]
    
    scell_markers = dplyr::filter(scell_markers, cluster %in% useful_cts)
    
    gene_markers = scell_markers$gene
    gene_markers_counts = table(gene_markers)
    
    #First get shared markers
    shared_markers = names(gene_markers_counts)[gene_markers_counts == length(useful_cts)]
    
    #Then unique markers
    scell_markers_unq = scell_markers %>%
      dplyr::filter(!gene %in% scell_markers$gene[duplicated(scell_markers$gene) == T]) %>%
      dplyr::slice(1:100) 
    
    unique_markers = scell_markers_unq %>% dplyr::select(gene) %>% pull()
    
    #Here define interesting spots: This must change
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    ct_scores = ct_scores_lt_all[[slide]]
    
    # Cells with less than .1 score for all cell-types associated with pattern
    pattern_ix_cols = grepl(ct_pattern, colnames(ct_scores),ignore.case = F)
    not_pattern_ids = rowSums(ct_scores[,pattern_ix_cols,drop=F] < .1) == length(useful_cts)
    not_pattern_ids = ct_scores$spot_id[not_pattern_ids]
    
    # Here we identify cells that have a chance of being the pattern cell
    ct_scores_mod = ct_scores %>%
      mutate(predicted_id = ifelse(grepl(ct_pattern,
                                         predicted_id,
                                         ignore.case = F) &
                                   prediction_score_max >= 0.1,
                                   predicted_id, "border"))
    
    # Here we take the border definition
    ct_scores_mod = ct_scores_mod %>% 
      mutate(predicted.id = ifelse(spot_id %in% not_pattern_ids,
                                   paste0("not_",ct_pattern), predicted_id))
    
    
    ct_scores_label = setNames(ct_scores_mod$predicted.id,ct_scores_mod$spot_id)
    
    visium_slide = AddMetaData(visium_slide,
                metadata = ct_scores_label[colnames(visium_slide)],col.name = "spot_label")
    
    Idents(visium_slide) = "spot_label"
    n_sample = length(unique(Idents(visium_slide)))
    set.seed(3)
    cells_analysed = SpatialDimPlot(visium_slide, label = TRUE, 
                   label.size = 0,stroke = 0,label.box = F) +
      scale_fill_manual(values = sample(RColorBrewer::brewer.pal(n = 8,"Accent"),n_sample))
    
    pdf(file = QC_out, width = 12, height = 12,onefile = TRUE)
    
    print(cells_analysed)
    
    # Slide coverage: This is to avoid dealing with sparsity in MISTy
    all_markers_init = c(shared_markers, unique_markers)
    
    slide_gex = visium_slide@assays$SCT@data
    
    slide_gex = as.matrix(slide_gex[all_markers_init[all_markers_init %in% rownames(slide_gex)],spot_ids])
    
    spot_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    # Coverage filter controlled by min proportion
    
    coverage_min = min(table(ct_scores_label[spot_ids])/ length(ct_scores_label[spot_ids]))
    coverage_min = coverage_min - (.5 * coverage_min)
    spot_coverage = names(spot_coverage[spot_coverage>=coverage_min])
    
    # Useful markers?
    shared_markers = shared_markers[shared_markers %in% spot_coverage]
    unique_markers = unique_markers[unique_markers %in% spot_coverage]
    
      #Merge them and make annotations
    all_markers = c(shared_markers, unique_markers)
    all_markers = all_markers[!grepl("-",all_markers)]
    all_markers_ann =  scell_markers %>%
      dplyr::filter(gene %in% all_markers) %>%
      mutate(shared = ifelse(gene %in% shared_markers,T,F)) %>%
      ungroup() %>% arrange(desc(shared),cluster) %>%
      dplyr::select(gene,cluster,shared)
    
    dotqc = Seurat::DotPlot(object = visium_slide,features = unique(all_markers_ann$gene),assay = "SCT") +
      coord_flip() + theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
    
    plot(dotqc)
    
    dev.off()
    
    #Which pathways?
    pathways = rownames(visium_slide@assays$progeny)
    
    #Which ligands?
    #To do
    ligands = rownames(visium_slide@assays$ligands)
    
    #Running MISTy
    
    #pathways
    ls = c(2,5,10,20,50,100)
    
    clear_cache()
    
    test_A1 = lapply(ls, run_mrkr_MISTy, 
                     seurat_visium_obj = visium_slide,
                     spot_ids = spot_ids,
                     intrinsic_markers = unique(all_markers),
                     ligands = NULL,
                     pathways = pathways,
                     out_dir_name = misty_out)
    
    #Getting optimal
    # get_optimal(out_dir_name = misty_out,ls = ls)
    # MISTy_out = MISTy_aggregator(results_folder = paste0(misty_out,"_optim"))
    # 
    # performance_file = paste0(misty_out,"_optim/performance.txt")
    # 
    # plot_MISTy_gene_path(all_markers_ann = all_markers_ann,
    #                      MISTy_out = MISTy_out,
    #                      pdf_out = pdf_out,
    #                      ligs = F,
    #                      performance = read_delim(performance_file, delim = " "))
    
    clear_cache()
    
    #ligands
    test_A2 = lapply(ls, run_mrkr_MISTy, 
                     seurat_visium_obj = visium_slide,
                     spot_ids = spot_ids,
                     intrinsic_markers = unique(all_markers),
                     ligands = ligands,
                     pathways = NULL,
                     out_dir_name = misty_out_ligs)
    
    #Getting optimal
    # get_optimal(out_dir_name = misty_out_ligs,ls = ls)
    # MISTy_out = MISTy_aggregator(results_folder = paste0(misty_out_ligs,"_optim"))
    # 
    # performance_file = paste0(misty_out_ligs,"_optim/performance.txt")
    # 
    # plot_MISTy_gene_path(all_markers_ann = all_markers_ann,
    #                      MISTy_out = MISTy_out,
    #                      pdf_out = pdf_out_ligs,
    #                      ligs = T,
    #                      performance = read_delim(performance_file, delim = " "))
    # 

  }
}



MISTy_marker_pipeline_down = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  ls = c(2,5,10,20,50,100)
  
  #Defining out files
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  misty_out = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_path_%s",
                      slide,slide,ct_pattern)
  
  misty_out_ligs = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_ligs_%s",
                           slide,slide,ct_pattern)
  
  pdf_out = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_path_%s.pdf",
                    slide,slide,ct_pattern)
  
  pdf_out_ligs = sprintf("./results/single_slide/%s/misty/%s_mrun_mrkr_ligs_%s.pdf",
                         slide,slide,ct_pattern)
  
  #Cell type scores from specific scRNA if data available
  scell_id = dplyr::filter(sample_dictionary, visium_ids == slide) %>%
    dplyr::select(scell_ids) %>% pull()
  
  visium_slide = readRDS(slide_file)
  
  if(scell_id %in% scell_marker_files){
    
    #Reading genesorter markers
    scell_markers = read.table(file = sprintf("./results/single_sample/top50/%s.top50.genes.txt",
                                              scell_id),
                               header = T, sep = "\t",
                               stringsAsFactors = F) %>%
      tidyr::pivot_longer(cols=colnames(.),"cluster",values_to = "gene") %>%
      dplyr::mutate(cluster = gsub("[.]","_",cluster))
    
    # Filter by most variable genes in the specific slide
    spec_vargns = all_var_genes %>% dplyr::filter(id == scell_id) %>%
      dplyr::select(value) %>% pull()
    
    scell_markers = dplyr::filter(scell_markers, gene %in% spec_vargns)
    
    # Continue processing of markers
    available_cts = unique(scell_markers$cluster)
    
    useful_cts = available_cts[grepl(pattern = ct_pattern, available_cts,
                                     ignore.case = F)]
    
    scell_markers = dplyr::filter(scell_markers, cluster %in% useful_cts)
    
    gene_markers = scell_markers$gene
    gene_markers_counts = table(gene_markers)
    
    #First get shared markers
    shared_markers = names(gene_markers_counts)[gene_markers_counts == length(useful_cts)]
    
    #Then unique markers
    scell_markers_unq = scell_markers %>%
      dplyr::filter(!gene %in% scell_markers$gene[duplicated(scell_markers$gene) == T]) %>%
      dplyr::slice(1:100) 
    
    unique_markers = scell_markers_unq %>% dplyr::select(gene) %>% pull()
    
    #Here define interesting spots: This must change
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    ct_scores = ct_scores_lt_all[[slide]]
    
    # Cells with less than .1 score for all cell-types associated with pattern
    pattern_ix_cols = grepl(ct_pattern, colnames(ct_scores),ignore.case = F)
    not_pattern_ids = rowSums(ct_scores[,pattern_ix_cols,drop=F] < .1) == length(useful_cts)
    not_pattern_ids = ct_scores$spot_id[not_pattern_ids]
    
    # Here we identify cells that have a chance of being the pattern cell
    ct_scores_mod = ct_scores %>%
      mutate(predicted_id = ifelse(grepl(ct_pattern,
                                         predicted_id,
                                         ignore.case = F) &
                                     prediction_score_max >= 0.1,
                                   predicted_id, "border"))
    
    # Here we take the border definition
    ct_scores_mod = ct_scores_mod %>% 
      mutate(predicted.id = ifelse(spot_id %in% not_pattern_ids,
                                   paste0("not_",ct_pattern), predicted_id))
    
    
    ct_scores_label = setNames(ct_scores_mod$predicted.id,ct_scores_mod$spot_id)
    
    # Slide coverage: This is to avoid dealing with sparsity in MISTy
    all_markers_init = c(shared_markers, unique_markers)
    
    slide_gex = visium_slide@assays$SCT@data
    
    slide_gex = as.matrix(slide_gex[all_markers_init[all_markers_init %in% rownames(slide_gex)],spot_ids])
    
    spot_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    # Coverage filter controlled by min proportion
    
    coverage_min = min(table(ct_scores_label[spot_ids])/ length(ct_scores_label[spot_ids]))
    coverage_min = coverage_min - (.5 * coverage_min)
    spot_coverage = names(spot_coverage[spot_coverage>=coverage_min])
    
    # Useful markers?
    shared_markers = shared_markers[shared_markers %in% spot_coverage]
    unique_markers = unique_markers[unique_markers %in% spot_coverage]
    
    #Merge them and make annotations
    all_markers = c(shared_markers, unique_markers)
    all_markers = all_markers[!grepl("-",all_markers)]
    all_markers_ann =  scell_markers %>%
      dplyr::filter(gene %in% all_markers) %>%
      mutate(shared = ifelse(gene %in% shared_markers,T,F)) %>%
      ungroup() %>% arrange(desc(shared),cluster) %>%
      dplyr::select(gene,cluster,shared)
    
    #pathways
    #Getting optimal
    get_optimal(out_dir_name = misty_out,ls = ls)
    MISTy_out = MISTy_aggregator(results_folder = paste0(misty_out,"_optim"))
    
    performance_file = paste0(misty_out,"_optim/performance.txt")
    
    plot_MISTy_gene_path(all_markers_ann = all_markers_ann,
                         MISTy_out = MISTy_out,
                         pdf_out = pdf_out,
                         ligs = F,
                         performance = read_delim(performance_file, delim = " "))
    
    
    #ligands
    #Getting optimal
    get_optimal(out_dir_name = misty_out_ligs,ls = ls)
    MISTy_out = MISTy_aggregator(results_folder = paste0(misty_out_ligs,"_optim"))
    
    performance_file = paste0(misty_out_ligs,"_optim/performance.txt")
    
    plot_MISTy_gene_path(all_markers_ann = all_markers_ann,
                         MISTy_out = MISTy_out,
                         pdf_out = pdf_out_ligs,
                         ligs = T,
                         performance = read_delim(performance_file, delim = " "))
    
  }
}



# CT-CT MISTy
MISTy_hpredictors_pipeline = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  
  #Defining out files
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  misty_out_ligs = sprintf("./results/single_slide/%s/misty/%s_hpredictors_%s",
                           slide,slide,ct_pattern)
  
  QC_out = sprintf("./results/single_slide/%s/misty/%s_%s_qc.pdf",
                   slide,slide,ct_pattern)
  
  all_markers_out = sprintf("./results/single_slide/%s/misty/%s_allmarkers_%s.rds",
                    slide,slide,ct_pattern)
  
  # General parameters
  hpredictors = c(readRDS("./markers/heart_predictors.rds"),c("NPPA","NPPB","CCR2"))
  visium_slide = readRDS(slide_file)
  
  ls = c(2,5,10,20,50,100)
  
  #Cell type scores from specific scRNA if data available
  scell_id = dplyr::filter(sample_dictionary, visium_ids == slide) %>%
    dplyr::select(scell_ids) %>% pull()
  
  if(scell_id %in% scell_marker_files){
    
    #Reading genesorter markers
    scell_markers = read.table(file = sprintf("./results/single_sample/top50/%s.top50.genes.txt",
                                              scell_id),
                               header = T, sep = "\t",
                               stringsAsFactors = F) %>%
      tidyr::pivot_longer(cols=colnames(.),"cluster",values_to = "gene") %>%
      dplyr::mutate(cluster = gsub("[.]","_",cluster))
    
    # Filter by most variable genes in the specific slide
    spec_vargns = all_var_genes %>% dplyr::filter(id == scell_id) %>%
      dplyr::select(value) %>% pull()
    
    scell_markers = dplyr::filter(scell_markers, gene %in% spec_vargns)
    
    # Continue processing of markers
    available_cts = unique(scell_markers$cluster)
    
    useful_cts = available_cts[grepl(pattern = ct_pattern, available_cts,
                                     ignore.case = F)]
    
    scell_markers = dplyr::filter(scell_markers, cluster %in% useful_cts)
    
    gene_markers = scell_markers$gene
    gene_markers_counts = table(gene_markers)
    
    #First get shared markers
    shared_markers = names(gene_markers_counts)[gene_markers_counts == length(useful_cts)]
    
    #Then unique markers
    scell_markers_unq = scell_markers %>%
      dplyr::filter(!gene %in% scell_markers$gene[duplicated(scell_markers$gene) == T]) %>%
      dplyr::slice(1:100) 
    
    unique_markers = scell_markers_unq %>% dplyr::select(gene) %>% pull()
    
    #Here define interesting spots: This must change
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    ct_scores = ct_scores_lt_all[[slide]]
    
    # Cells with less than .1 score for all cell-types associated with pattern
    pattern_ix_cols = grepl(ct_pattern, colnames(ct_scores),ignore.case = F)
    not_pattern_ids = rowSums(ct_scores[,pattern_ix_cols,drop=F] < .1) == length(useful_cts)
    not_pattern_ids = ct_scores$spot_id[not_pattern_ids]
    
    # Here we identify cells that have a chance of being the pattern cell
    ct_scores_mod = ct_scores %>%
      mutate(predicted_id = ifelse(grepl(ct_pattern,
                                         predicted_id,
                                         ignore.case = F) &
                                     prediction_score_max >= 0.1,
                                   predicted_id, "border"))
    
    # Here we take the border definition
    ct_scores_mod = ct_scores_mod %>% 
      mutate(predicted.id = ifelse(spot_id %in% not_pattern_ids,
                                   paste0("not_",ct_pattern), predicted_id))
    
    
    ct_scores_label = setNames(ct_scores_mod$predicted.id,ct_scores_mod$spot_id)
    
    visium_slide = AddMetaData(visium_slide,
                               metadata = ct_scores_label[colnames(visium_slide)],col.name = "spot_label")
    
    Idents(visium_slide) = "spot_label"
    n_sample = length(unique(Idents(visium_slide)))
    set.seed(3)
    cells_analysed = SpatialDimPlot(visium_slide, label = TRUE, 
                                    label.size = 0,stroke = 0,label.box = F) +
      scale_fill_manual(values = sample(RColorBrewer::brewer.pal(n = 8,"Accent"),n_sample))
    
    pdf(file = QC_out, width = 12, height = 12,onefile = TRUE)
    
    print(cells_analysed)
    
    # Slide coverage: This is to avoid dealing with sparsity in MISTy
    all_markers_init = c(shared_markers, unique_markers)
    
    slide_gex = visium_slide@assays$SCT@data
    
    slide_gex = as.matrix(slide_gex[all_markers_init[all_markers_init %in% rownames(slide_gex)],spot_ids])
    
    spot_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    # Coverage filter controlled by min proportion
    
    coverage_min = min(table(ct_scores_label[spot_ids])/ length(ct_scores_label[spot_ids]))
    coverage_min = coverage_min - (.5 * coverage_min)
    spot_coverage = names(spot_coverage[spot_coverage>=coverage_min])
    
    # Useful markers?
    shared_markers = shared_markers[shared_markers %in% spot_coverage]
    unique_markers = unique_markers[unique_markers %in% spot_coverage]
    
    #Merge them and make annotations
    all_markers = c(shared_markers, unique_markers)
    all_markers = all_markers[!grepl("-",all_markers)]
    all_markers_ann =  scell_markers %>%
      dplyr::filter(gene %in% all_markers) %>%
      mutate(shared = ifelse(gene %in% shared_markers,T,F)) %>%
      ungroup() %>% arrange(desc(shared),cluster) %>%
      dplyr::select(gene,cluster,shared)
    
    saveRDS(all_markers_ann,file = all_markers_out)
    
    dotqc = Seurat::DotPlot(object = visium_slide,features = unique(all_markers_ann$gene),assay = "SCT") +
      coord_flip() + theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
    
    plot(dotqc)
    
    dev.off()
    
    #Which ligands?
    slide_gex = visium_slide@assays$SCT@data
    slide_gex = as.matrix(slide_gex[hpredictors[hpredictors %in% rownames(slide_gex)],])
    predictor_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    ligands = names(predictor_coverage[predictor_coverage>0.1])
    
    #Running MISTy
    
    #pathways
    ls = c(2,5,10,20,50,100)
    
    clear_cache()
    
    test_A1 = lapply(ls, run_mrkr_MISTy, 
                     seurat_visium_obj = visium_slide,
                     spot_ids = spot_ids,
                     intrinsic_markers = unique(all_markers),
                     ligands = unique(ligands),
                     pathways = NULL,
                     out_dir_name = misty_out_ligs)
    
  }
}

#
#
#
MISTy_marker_pipeline_down_v2 = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  
  #Defining out files
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  all_markers_out = sprintf("./results/single_slide/%s/misty/%s_allmarkers_%s.rds",
                            slide,slide,ct_pattern)
  
  misty_out = sprintf("./results/single_slide/%s/misty/%s_hpredictors_%s",
                      slide,slide,ct_pattern)
  
  pdf_out = sprintf("./results/single_slide/%s/misty/%s_hpredictors_summary_%s.pdf",
                    slide,slide,ct_pattern)
  
  
  visium_slide = readRDS(slide_file)
  
  all_markers_ann = readRDS(all_markers_out)
  
  #Genes to perform differential expression analysis
  predicted_markers = unique(all_markers_ann$gene)
  
  #Prepare visium object
  ct_scores = ct_scores_lt_all[[slide]]
  
  # Cells with less than .1 score for all cell-types associated with pattern
  pattern_ix_cols = grepl(ct_pattern, colnames(ct_scores),ignore.case = F)
  
  useful_cts = gsub("prediction_score_","",colnames(ct_scores)[pattern_ix_cols])
  not_pattern_ids = rowSums(ct_scores[,pattern_ix_cols,drop=F] < .1) == length(useful_cts)
  not_pattern_ids = ct_scores$spot_id[not_pattern_ids]
  
  # Here we identify cells that have a chance of being the pattern cell
  ct_scores_mod = ct_scores %>%
    mutate(predicted_id = ifelse(grepl(ct_pattern,
                                       predicted_id,
                                       ignore.case = F) &
                                   prediction_score_max >= 0.1,
                                 predicted_id, "border"))
  
  # Here we take the border definition
  ct_scores_mod = ct_scores_mod %>% 
    mutate(predicted.id = ifelse(spot_id %in% not_pattern_ids,
                                 paste0("not_",ct_pattern), predicted_id))
  
  
  ct_scores_label = setNames(ct_scores_mod$predicted.id,ct_scores_mod$spot_id)
  
  visium_slide = AddMetaData(visium_slide,
                             metadata = ct_scores_label[colnames(visium_slide)],col.name = "spot_label")
  
  Idents(visium_slide) = "spot_label"
  
  # New annotation
  spec_markrs = FindAllMarkers(object = visium_slide,test.use = "wilcox",
                               features = predicted_markers,logfc.threshold = 0.05,only.pos = T)
  
  spec_markrs = spec_markrs %>% dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(cluster,-avg_logFC) %>%
    dplyr::mutate(ct = as.character(cluster)) %>%
    dplyr::select(ct, gene)
  
  #shared_expression = table(spec_markrs$gene)
  
  #shared_genes = names(shared_expression[shared_expression>1])
  
  # spec_markrs = spec_markrs %>%
  #   dplyr::mutate(ct = as.character(cluster)) %>%
  #   dplyr::mutate(ct = ifelse(gene %in% shared_genes,
  #                             "shared",ct)) %>%
  #   dplyr::select(ct,gene) %>% unique()
  
  # Getting optim results
  
  MISTy_out = MISTy_aggregator(results_folder = paste0(misty_out,"_optim"))
  
  # Plotting essential
  
  # First colors and order
  set.seed(117)
  
  color_df = unique(spec_markrs$ct)
  color_v = sample((RColorBrewer::brewer.pal(n = 8,
                                             "Dark2")), length(color_df))
  
  color_df = tibble(ct = color_df, ccolor = color_v)
  
  shared_expression = table(spec_markrs$gene)
  
  shared_genes = names(shared_expression[shared_expression>1])
  
  color_df = left_join(spec_markrs, color_df) %>% arrange(ct)
  
  marker_plt = ggplot(color_df, aes(x = factor(gene,
                                               levels = unique(gene)),
                                    y = ct,
                                    fill = ct)) +
    geom_tile() + theme_minimal() + 
    theme(
      axis.title = element_text(size=11),
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5))
  
  # Overall performance
  performance = read_delim(paste0(misty_out,"_optim/performance.txt"), 
                           delim = " ")
  
  performance_plt = performance %>% tidyr::pivot_longer(cols = -target) %>%
    dplyr::filter(grepl("R2",name)) %>%
    dplyr::filter(target %in% color_df$gene) %>%
    dplyr::mutate(target = factor(target,
                                  levels = unique(color_df$gene))) %>%
    ggplot(aes(fill = name, y = target, x = value)) +
    geom_bar(stat = "identity",position="dodge") + 
    theme_minimal() #+
    #theme(axis.text.y = element_text(color = color_df$ccolor))
  
  R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
    mutate(target = targets) %>%
    dplyr::filter(target %in% color_df$gene) %>%
    pivot_longer(cols = -target, 
                 names_to = "name", 
                 values_to = "value") %>% arrange(desc(value)) %>%
    mutate(target = factor(target,
                           levels = unique(color_df$gene)))
  
  impr_plot = ggplot(R2_impr) +
    geom_point(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Change in variance explained") +
    xlab("Target") +
    #ylim(c(-5, 25)) +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  #Look at importances of views
  contribution_df = MISTy_out$coefs %>% dplyr::filter(target %in% R2_impr$target) %>%
    dplyr::filter(target %in% color_df$gene)
  
  coefs_plot = ggplot(contribution_df) + 
    geom_col(aes(x=factor(target, levels = unique(color_df$gene)), 
                 y=value, group=view, fill=view)) +
    xlab("Target") +
    ylab("Contribution") +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          legend.text = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     vjust = 0.5)) +
    scale_fill_brewer(palette="Dark2",
                      labels = c("Intrinsic","Para pathway")) +
    labs(fill = "View")
  
  
  #Intra
  importance_intra = tidyr::gather(MISTy_out$importance[[1]], "Predicted",
                                   "Importance", -Predictor) %>%
    left_join(all_markers_ann, by = c("Predicted"="gene")) %>%
    mutate(PredictedCell = cluster) %>% 
    select(Predictor,Predicted,
           Importance,PredictedCell) %>%
    left_join(all_markers_ann, by = c("Predictor"="gene")) %>%
    mutate(PredictorCell = cluster) %>% 
    select(Predictor,Predicted,Importance,
           PredictorCell,
           PredictedCell) %>%
    arrange(PredictedCell,PredictorCell,Predicted,Predictor) %>%
    dplyr::filter(Predicted %in% color_df$gene)
  
  importance_intra = importance_intra %>%
    mutate(Predicted = factor(Predicted,
                              levels = unique(color_df$gene)))
  
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
  
  #Para
  
  importance_para = tidyr::gather(MISTy_out$importance[[2]], "Predicted",
                                  "Importance", -Predictor) %>%
    dplyr::filter(Predicted %in% color_df$gene) %>%
    arrange(-Importance) %>%
    left_join(color_df,by = c("Predicted" = "gene"))
  
  # Identify top n predictors of all markers
  
  best_predictors = importance_para %>% group_by(ct,Predictor) %>%
    summarise(mean_importance = mean(Importance)) %>%
    ungroup() %>%
    arrange(ct,-mean_importance) %>%
    dplyr::filter(mean_importance>1) %>%
    unique()
  
  bp_plt = ggplot(best_predictors, aes(x = factor(Predictor,
                                                  levels = unique(Predictor)),
                                       y = ct,
                                       fill = mean_importance)) +
    geom_tile() + theme_minimal() + 
    theme(
      axis.title = element_text(size=11),
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5))
  
  importance_para = importance_para %>%
    dplyr::filter(Predictor %in%
                    best_predictors$Predictor) %>%
    unique() %>%
    dplyr::select(Predictor, Predicted, Importance) %>%
    dplyr::mutate(Predictor = factor(Predictor,
                                     levels = unique(best_predictors$Predictor)),
                  Predicted = factor(Predicted,
                                     levels = unique(color_df$gene)))
  
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
    xlab("Para ligands")
  
  
  pdf(file = pdf_out, width = 15, height = 9,onefile = TRUE)
  plot(marker_plt)
  plot(performance_plt)
  plot(impr_plot)
  plot(coefs_plot)
  plot(importance_intra_plt)
  plot(bp_plt)
  plot(importance_para_plt)
  dev.off()
  
  
  out_res = list("marker_dictionary" = spec_markrs,
                 "intra_importance" = unique(importance_intra[,1:3]),
                 "para_importance" = unique(importance_para),
                 "para_summary" = best_predictors
                 )
  
  return(out_res)
}

# CT-CT MISTy
MISTy_cytokines_pipeline = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  
  #Defining out files
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  misty_out_ligs = sprintf("./results/single_slide/%s/misty/%s_cytok_%s",
                           slide,slide,ct_pattern)
  
  QC_out = sprintf("./results/single_slide/%s/misty/%s_%s_cytokqc.pdf",
                   slide,slide,ct_pattern)
  
  all_markers_out = sprintf("./results/single_slide/%s/misty/%s_allmarkers_cytok_%s.rds",
                            slide,slide,ct_pattern)
  
  # General parameters
  gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")
  cytokines = gene_sets$MSIGDB_KEGG$KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
  hpredictors = cytokines
  visium_slide = readRDS(slide_file)
  
  ls = c(2,5,10,20,50,100)
  
  #Cell type scores from specific scRNA if data available
  scell_id = dplyr::filter(sample_dictionary, visium_ids == slide) %>%
    dplyr::select(scell_ids) %>% pull()
  
  if(scell_id %in% scell_marker_files){
    
    #Reading genesorter markers
    scell_markers = read.table(file = sprintf("./results/single_sample/top50/%s.top50.genes.txt",
                                              scell_id),
                               header = T, sep = "\t",
                               stringsAsFactors = F) %>%
      tidyr::pivot_longer(cols=colnames(.),"cluster",values_to = "gene") %>%
      dplyr::mutate(cluster = gsub("[.]","_",cluster))
    
    # Filter by most variable genes in the specific slide
    spec_vargns = all_var_genes %>% dplyr::filter(id == scell_id) %>%
      dplyr::select(value) %>% pull()
    
    scell_markers = dplyr::filter(scell_markers, gene %in% spec_vargns)
    
    # Continue processing of markers
    available_cts = unique(scell_markers$cluster)
    
    useful_cts = available_cts[grepl(pattern = ct_pattern, available_cts,
                                     ignore.case = F)]
    
    scell_markers = dplyr::filter(scell_markers, cluster %in% useful_cts)
    
    gene_markers = scell_markers$gene
    gene_markers_counts = table(gene_markers)
    
    #First get shared markers
    shared_markers = names(gene_markers_counts)[gene_markers_counts == length(useful_cts)]
    
    #Then unique markers
    scell_markers_unq = scell_markers %>%
      dplyr::filter(!gene %in% scell_markers$gene[duplicated(scell_markers$gene) == T]) %>%
      dplyr::slice(1:100) 
    
    unique_markers = scell_markers_unq %>% dplyr::select(gene) %>% pull()
    
    #Here define interesting spots: This must change
    
    spot_ids = get_CTcoordinates(ct_pattern = ct_pattern,
                                 slide = slide, 
                                 max_filter = 0.1)
    
    ct_scores = ct_scores_lt_all[[slide]]
    
    # Cells with less than .1 score for all cell-types associated with pattern
    pattern_ix_cols = grepl(ct_pattern, colnames(ct_scores),ignore.case = F)
    not_pattern_ids = rowSums(ct_scores[,pattern_ix_cols,drop=F] < .1) == length(useful_cts)
    not_pattern_ids = ct_scores$spot_id[not_pattern_ids]
    
    # Here we identify cells that have a chance of being the pattern cell
    ct_scores_mod = ct_scores %>%
      mutate(predicted_id = ifelse(grepl(ct_pattern,
                                         predicted_id,
                                         ignore.case = F) &
                                     prediction_score_max >= 0.1,
                                   predicted_id, "border"))
    
    # Here we take the border definition
    ct_scores_mod = ct_scores_mod %>% 
      mutate(predicted.id = ifelse(spot_id %in% not_pattern_ids,
                                   paste0("not_",ct_pattern), predicted_id))
    
    
    ct_scores_label = setNames(ct_scores_mod$predicted.id,ct_scores_mod$spot_id)
    
    visium_slide = AddMetaData(visium_slide,
                               metadata = ct_scores_label[colnames(visium_slide)],col.name = "spot_label")
    
    Idents(visium_slide) = "spot_label"
    n_sample = length(unique(Idents(visium_slide)))
    set.seed(3)
    cells_analysed = SpatialDimPlot(visium_slide, label = TRUE, 
                                    label.size = 0,stroke = 0,label.box = F) +
      scale_fill_manual(values = sample(RColorBrewer::brewer.pal(n = 8,"Accent"),n_sample))
    
    pdf(file = QC_out, width = 12, height = 12,onefile = TRUE)
    
    print(cells_analysed)
    
    # Slide coverage: This is to avoid dealing with sparsity in MISTy
    all_markers_init = c(shared_markers, unique_markers)
    
    slide_gex = visium_slide@assays$SCT@data
    
    slide_gex = as.matrix(slide_gex[all_markers_init[all_markers_init %in% rownames(slide_gex)],spot_ids])
    
    spot_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    # Coverage filter controlled by min proportion
    
    coverage_min = min(table(ct_scores_label[spot_ids])/ length(ct_scores_label[spot_ids]))
    coverage_min = coverage_min - (.5 * coverage_min)
    spot_coverage = names(spot_coverage[spot_coverage>=coverage_min])
    
    # Useful markers?
    shared_markers = shared_markers[shared_markers %in% spot_coverage]
    unique_markers = unique_markers[unique_markers %in% spot_coverage]
    
    #Merge them and make annotations
    all_markers = c(shared_markers, unique_markers)
    all_markers = all_markers[!grepl("-",all_markers)]
    all_markers_ann =  scell_markers %>%
      dplyr::filter(gene %in% all_markers) %>%
      mutate(shared = ifelse(gene %in% shared_markers,T,F)) %>%
      ungroup() %>% arrange(desc(shared),cluster) %>%
      dplyr::select(gene,cluster,shared)
    
    saveRDS(all_markers_ann,file = all_markers_out)
    
    dotqc = Seurat::DotPlot(object = visium_slide,features = unique(all_markers_ann$gene),assay = "SCT") +
      coord_flip() + theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
    
    plot(dotqc)
    
    dev.off()
    
    #Which ligands?
    slide_gex = visium_slide@assays$SCT@data
    slide_gex = as.matrix(slide_gex[hpredictors[hpredictors %in% rownames(slide_gex)],])
    predictor_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    
    ligands = names(predictor_coverage[predictor_coverage>0.01])
    
    #Running MISTy
    
    #pathways
    ls = c(2,5,10,20,50,100)
    
    clear_cache()
    
    test_A1 = lapply(ls, run_mrkr_MISTy, 
                     seurat_visium_obj = visium_slide,
                     spot_ids = spot_ids,
                     intrinsic_markers = unique(all_markers),
                     ligands = unique(ligands),
                     pathways = NULL,
                     out_dir_name = misty_out_ligs)
    
  }
}