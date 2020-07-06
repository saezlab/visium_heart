# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' MISTy runs in fibrotic area and borderzone
#' Para-views are pathways, ECM and cytokines
#' 
#' 
#' 
#' 

library(tidyverse)
library(Seurat)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

# GLOBAL VARIABLES

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

# CT scores coordinates
get_CTcoordinates = function(ct_pattern, slide, max_filter = 0){
  ct_scores = ct_scores_lt_all[[slide]]
  ixs = grepl(ct_pattern,ct_scores$predicted_id,ignore.case = F)
  ct_scores = ct_scores[ixs,]
  ixs = ct_scores$prediction_score_max >= max_filter
  ct_ids = ct_scores$spot_id[ixs]
  return(ct_ids)
}

# MISTY MAIN PIPELINE

MISTy_heart_pipeline = function(slide, ct_pattern){
  
  print(slide)
  print(ct_pattern)
  
  #Defining out files
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  QC_out = sprintf("./results/single_slide/%s/misty/%s_%s_qc.pdf",
                   slide,slide,ct_pattern)
  
  misty_ecm_out = sprintf("./results/single_slide/%s/misty/%s_%s_ecm",
                          slide,slide,ct_pattern)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_%s_ctk",
                            slide,slide,ct_pattern)
  
  misty_path_out = sprintf("./results/single_slide/%s/misty/%s_%s_path",
                           slide,slide,ct_pattern)
  
  all_markers_out = sprintf("./results/single_slide/%s/misty/%s_allmarkers_%s.rds",
                            slide,slide,ct_pattern)
  
  # Reading the slide
  visium_slide = readRDS(slide_file)
  DefaultAssay(visium_slide) = "SCT"
  
  
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
    
    #######################
    #Running MISTy
    
    ls = c(2,5,10,20,50,100)
    
    # This is PROGENy high level analysis
    
    test_path = lapply(ls,para_ppln_seurat,
                       visium_slide = visium_slide,
                       intra_assay = "SCT",
                       intra_features = unique(all_markers),
                       para_assay = "progeny",
                       para_features = NULL,
                       spot_ids = NULL,
                       out_alias = misty_path_out)
    
    # This is ECM analysis
    
    # Prepare MISTy parameters
    hpredictors = c(readRDS("./markers/heart_predictors.rds"),
                    c("NPPA","NPPB","CCR2"))
    #Which ECM prots?
    slide_gex = visium_slide@assays$SCT@data
    slide_gex = as.matrix(slide_gex[hpredictors[hpredictors %in% rownames(slide_gex)],])
    predictor_coverage = rowSums(slide_gex > 0)/ncol(slide_gex)
    ligands = names(predictor_coverage[predictor_coverage>0.1])
    
    test_ecm = lapply(ls,para_ppln_seurat,visium_slide = visium_slide,
                       intra_assay = "SCT",
                       intra_features = unique(all_markers),
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
                        intra_features = unique(all_markers),
                        para_assay = "SCT",
                        para_features = ligands,
                        spot_ids = NULL,
                        out_alias = misty_cytok_out)
    
    return(NULL)
  }
}


# Generate optims
wrap_optim = function(slide, ct_pattern){
  print(slide)
  print(ct_pattern)
  
  misty_ecm_out = sprintf("./results/single_slide/%s/misty/%s_%s_ecm",
                          slide,slide,ct_pattern)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_%s_ctk",
                            slide,slide,ct_pattern)
  
  misty_path_out = sprintf("./results/single_slide/%s/misty/%s_%s_path",
                           slide,slide,ct_pattern)
  
  ls = c(2,5,10,20,50,100)
  a = sapply(c(misty_ecm_out,misty_cytok_out,misty_path_out),
             get_optimal,ls = ls)
}

## MAIN CODE ##

selected_tests = readRDS(file = "./results/single_slide/MISTy_selection.rds") %>%
  dplyr::filter(slide %in% c("157772","157781"))

run_out = selected_tests[,1:2] %>%
  pmap(.,MISTy_heart_pipeline)

optim = selected_tests[,1:2] %>%
  pmap(.,wrap_optim)

wrap_optim(slide = "157781",ct_pattern = "Cardiomyocytes")
