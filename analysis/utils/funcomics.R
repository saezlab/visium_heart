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
library(tidyverse)
library(Seurat)
library(dorothea)
library(viper)
library(progeny)
library(decoupleR)


#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param confidence_lbls: confidence of regulons picked
#' @return A Seurat object with TF activities in dorothea assay
add_tf_activities <- function(visium_slide, 
                              species = "human",
                              confidence_lbls = c("A","B","C","D"),
                              verbose = FALSE,
                              assay = "RNA"){
  if (species == "mouse"){
    
    data(dorothea_mm, package = "dorothea")
    regulons <- dorothea_mm
    
  } else if (species == "human"){
    
    data(dorothea_hs, package = "dorothea")
    regulons <- dorothea_hs
  }
  
  regulons = regulons %>% 
    filter(confidence %in% confidence_lbls)
  
  tf_act_mat <- run_viper(input = as.matrix(GetAssayData(visium_slide, assay = assay)),
                          regulons = regulons, 
                          options = list(nes = TRUE, 
                                         method = "scale", minsize = 4, 
                                         eset.filter = FALSE,
                                         verbose = verbose))
  
  ## Repeated regulon RFXAP
  # tf_act_mat = tf_act_mat[!duplicated(tf_act_mat),]
  visium_slide[['dorothea']] = CreateAssayObject(data = tf_act_mat)
  
  return(visium_slide)
}




#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param top: number of genes used in footprint
#' @return A Seurat object with pathway activities in progeny assay
add_path_activities <- function(visium_slide, 
                               species = "human",
                               top = 500, 
                               verbose = F,
                               assay = "RNA"){
  
  if(species == "mouse"){
    model <- progeny::getModel(organism = "Mouse", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- scale(t(progeny_scores))
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
    #progeny_scores <- progeny::progeny(expr = GetAssayData(visium_slide, assay = assay),
    #                                  scale=TRUE, 
    #                                  organism="Mouse", 
    #                                  top=top, 
    #                                  perm=1, 
    #                                  verbose = verbose)
    
    #visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
  }else if(species == "human"){
    
    model <- progeny::getModel(organism = "Human", top = top)
    common_genes <- intersect(rownames(GetAssayData(visium_slide, assay = assay)), rownames(model))
    progeny_scores <- t(model)[, common_genes] %*% GetAssayData(visium_slide, assay = assay)[common_genes, ]
    progeny_scores <- scale(t(progeny_scores))
    
    visium_slide[['progeny']] <- CreateAssayObject(counts = t(progeny_scores))
    
    #progeny_scores <- progeny::progeny(expr = GetAssayData(visium_slide, assay = assay),
    #                         scale=TRUE, 
    #                         organism="Human", 
    #                         top=top, 
    #                         perm=1,
    #                         verbose = verbose)
    
    #visium_slide[['progeny']] <- CreateAssayObject(data = t(progeny_scores))
    
  }
  
  return(visium_slide)
}

# Get module scores function -----------------------------------------------------------------------------------
#' Generates a module score matrix
#' 
#' @param visium_slide: Seurat object
#' @param MS_regulon: list of gene sets
#' @return a matrix with module scores ready to be added as an assay
getTF_matrix_MS <- function(visium_slide, 
                            MS_regulon,
                            assay = "RNA",
                            module_name = "user_gsets") {
  
  names_vect <- gsub("[.]","_",names(MS_regulon))
  names_vect <- gsub("-","_",names_vect)
  
  tf_act_mat <- AddModuleScore(visium_slide,
                               assay = assay,
                               features = MS_regulon,
                               name = paste0(names_vect,"__"))
  
  tf_act_mat <- tf_act_mat@meta.data
  
  cell_ids <- rownames(tf_act_mat)
  calculated_regulons <- colnames(tf_act_mat)[grepl("__", colnames(tf_act_mat))]
  
  tf_act_mat <- tf_act_mat[,calculated_regulons]
  
  colnames(tf_act_mat) <- unlist(map(strsplit(colnames(tf_act_mat),split = "__"), function(x) x[1]))
  
  rownames(tf_act_mat) <- cell_ids
  
  tf_act_mat <- t(as.matrix(tf_act_mat))
  
  visium_slide[[module_name]] <- CreateAssayObject(data = tf_act_mat)

  return(visium_slide)
}

#' @param visium_slide: Seurat object with SCT assay
#' @param species: human or mouse, will extract regulons based on this
#' @param confidence_lbls: confidence of regulons picked for dorothea
#' @param top: number of genes used in footprint
#' @return A Seurat object with pathway activities in progeny assay
add_funcomics <- function(visium_slide,
                         species = "human",
                         confidence_lbls = c("A","B","C","D"),
                         top = 500,
                         marker_dictionary = NULL,
                         module_name = NULL,
                         verbose = FALSE,
                         assay = "RNA"){
  
  visium_slide <- add_tf_activities(visium_slide = visium_slide,
                                   species = species,
                                   confidence_lbls = confidence_lbls,
                                   verbose = verbose,
                                   assay = assay)
  
  visium_slide <- add_path_activities(visium_slide = visium_slide,
                                     species = species,
                                     top = top,
                                     verbose = verbose,
                                     assay = assay)
  
  if(!is.null(marker_dictionary)){
    visium_slide <- getTF_matrix_MS(visium_slide = visium_slide,
                                    MS_regulon = marker_dictionary,
                                    assay = assay,
                                    module_name = module_name)
  }
  
  return(visium_slide)
}

# Get module scores function -----------------------------------------------------------------------------------
#' Generates a module score matrix
#' 
#' @param visium_slide: Seurat object
#' @param MS_regulon: list of gene sets
#' @return a visium object
get_wmean_score <- function(visium_slide, 
                            network,
                            assay = "RNA",
                            module_name = "user_gsets") {

  
  gset_mat <- run_wmean(network = network,
                        mat = GetAssayData(visium_slide, assay = assay),
                        .source = "source",
                        .target = "target",
                        .mor = "mor",
                        .likelihood = "likelihood",
                        times = 100)

  gset_mat <- gset_mat %>% 
    dplyr::filter(statistic == "norm_wmean") %>%
    dplyr::select(-c("statistic","p_value")) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()
  
  visium_slide[[module_name]] <- CreateAssayObject(data = gset_mat)
  
  return(visium_slide)
}

