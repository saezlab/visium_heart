# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Runs funcomic tools as selected by user
#' 1) TF activities (dorothea)
#' 2) Pathway activities (progeny)
#' 3) Cell_type scores (if a data frame of markers/stats provided)
#' 4) Functional scores (if a data frame of gene sets is provided)
#' 
#' It must work regardless the species
#' 

#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param confidence_lbls: confidence of regulons picked
#' @return A Seurat object with TF activities in dorothea assay
add_tf_activities = function(visium_slide, 
                             species = "human",
                             confidence_lbls = c("A","B","C","D"),
                             verbose = FALSE){
  require(Seurat)
  require(dorothea)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(dplyr)
  require(viper)
  
  if (species == "mouse"){
    
  data(dorothea_mm, package = "dorothea")
  regulons = dorothea_mm
    
  } else if (species == "human"){
    
    data(dorothea_hs, package = "dorothea")
    regulons = dorothea_hs
  }
  
  regulons = regulons %>% 
    filter(confidence %in% confidence_lbls)
  
  tf_act_mat = run_viper(input = as.matrix(visium_slide[["SCT"]]@data), 
                         regulons = regulons, 
                         options = list(nes = TRUE, 
                                        method = "scale", minsize = 4, 
                                        eset.filter = FALSE,
                                        verbose = verbose)
                         )
  
  ## Repeated regulon RFXAP
  # tf_act_mat = tf_act_mat[!duplicated(tf_act_mat),]
  visium_slide[['dorothea']] = CreateAssayObject(data = tf_act_mat)
  
  return(visium_slide)
}


#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param top: number of genes used in footprint
#' @return A Seurat object with pathway activities in progeny assay
add_path_activities = function(visium_slide, 
                             species = "human",
                             top = 500, verbose = F){
  require(Seurat)
  require(progeny)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(dplyr)
  
  if(species == "mouse"){
    
    progeny_scores = progeny::progeny(expr = as.matrix(visium_slide[["SCT"]]@data),
                                      scale=TRUE, 
                                      organism="Mouse", top=top, perm=1, verbose = verbose)
    
    visium_slide[['progeny']] = CreateAssayObject(counts = t(progeny_scores))
    
  }else if(species == "human"){
    
    progeny_scores = progeny(expr = as.matrix(visium_slide[["SCT"]]@data),
                             scale=TRUE, 
                             organism="Human", top=top, perm=1)
    
    visium_slide[['progeny']] = CreateAssayObject(data = t(progeny_scores))
    
  }
  
  return(visium_slide)
}

#' @param visium_slide: Seurat object with SCT assay
#' @param marker_df: data frame w 3 columns (gene,size_effect,cell_type), previously
#' filtered marker genes
#' @return A Seurat object with cell type scores in ctscores assay
add_ctscores = function(visium_slide, 
                        marker_df){
  require(Seurat)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(dplyr)
  
  colnames(marker_df) = c("gene","size_effect","ct")
  marker_scores = pivot_wider(marker_df,id_cols = gene,
                              values_from = size_effect,names_from = ct)
  
  marker_scores_mat = as.matrix(marker_scores[,-1])
  rownames(marker_scores_mat) = marker_scores[[1]]
  
  colnames(marker_scores_mat) = gsub(" ","_",colnames(marker_scores_mat)) 
  
  #Linear combination
  exprdata = as.matrix(visium_slide[["SCT"]]@data)
  genes_in_markers = rownames(exprdata) %in% rownames(marker_scores_mat)
  exprdata = exprdata[genes_in_markers,]
  marker_scores_mat = marker_scores_mat[rownames(exprdata),]
  
  ct_scores = scale(t(t(marker_scores_mat) %*% exprdata))
  
  visium_slide[['ctscores']] = CreateAssayObject(counts = t(ct_scores))
  
  return(visium_slide)
}



#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param confidence_lbls: confidence of regulons picked for dorothea
#' @param top: number of genes used in footprint
#' @return A Seurat object with pathway activities in progeny assay
add_funcomics = function(visium_slide,
                         species = "human",
                         confidence_lbls = c("A","B","C","D"),
                         top = 500,
                         marker_df = NULL,
                         verbose = FALSE){
  
  visium_slide = add_tf_activities(visium_slide = visium_slide,
                                   species = species,
                                   confidence_lbls = confidence_lbls,
                                   verbose = verbose)
  
  visium_slide = add_path_activities(visium_slide = visium_slide,
                                     species = species,
                                     top = top,
                                     verbose = verbose)
  
  if(!is.null(marker_df)){
    visium_slide = add_ctscores(visium_slide = visium_slide,
                                marker_df = marker_df)
  }
  
  return(visium_slide)
}
                         



