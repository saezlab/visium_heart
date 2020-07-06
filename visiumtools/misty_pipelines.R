# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Defines general misty pipelines
#' 
#' Intra - para in all slide : Assuming that data is inside a slot 
#' (will generalize first to matrices)
#' 
#' Intra - para in specific spots 
#' (fitting in everything and then only getting the output for selected spots):
#' 
#' 
require(Seurat)
require(readr)
require(tidyr)
require(dplyr)
require(purrr)
require(tibble)
require(ggplot2)
require(stringr)
require(cowplot)
library(MISTy)
require(future)

#' Runs MISTy classic paraview pipeline 
#' 
#' Warning: Feature IDs can't have "-" symbol
#' Warning: rownames must be the same in intra_df, para_df and geometry
#'
#' @param intra_df: feature data_frame with IDs as rows and features as columns
#' @param intra_features: features of the intraview to be used, if NULL uses all columns
#' @param para_df: feature data_frame with IDs as rows and features as columns
#' @param para_features: features of the paraview to be used, if NULL uses all columns
#' @param geometry: a data frame with IDs as rows two coordinates as columns 
#' @param l: radius parameter
#' @param spot_ids: spot IDs to fit MISTy if null all rows are used
#' @param out_alias: folder name to be used in all MISTy outputs
para_pipeline = function(intra_df, 
                               intra_features = NULL,
                               para_df,
                               para_features = NULL,
                               geometry,
                               l,
                               spot_ids = NULL,
                               out_alias = "default"){
  
  plan(multiprocess, workers = 4)
  
  clear_cache()
  
  if(is.null(intra_features)){
    intra_features = colnames(intra_df)
  }
  
  if(is.null(para_features)){
    para_features = colnames(para_df)
  }
  
  if(is.null(spot_ids)){
    spot_ids = rownames(intra_df)
  }
  
  # Defining useful data intra
  intra_df = intra_df[spot_ids,intra_features]
  colnames(intra_df) = gsub("-","_", colnames(intra_df))
  
  views_main = create_initial_view(intra_df, 
                                   unique.id = "intra")
  
  
  # Defining useful data para
  para_df = para_df[rownames(geometry),para_features]
  
  colnames(para_df) = gsub("-","_", colnames(para_df))
  
  views_para = create_initial_view(para_df, 
                                   unique.id = paste0("para_",l^2)) %>% 
    add_paraview(geometry, l^2)
  
  # Fetching actual para info to be used
  
  # Spot specific view comes from the view above
  data_red = views_para[[3]]$data
  rownames(data_red) = rownames(para_df) #we named rows just for easy access
  data_red = data_red[rownames(intra_df),]
  
  # Define frankenstein view
  views = views_main %>% 
    add_views(create_view(paste0("para_",l^2),
                          data_red))
  
  MISTy_run = run_misty(views,paste0(out_alias,"_",l^2))
  
}


para_ppln_seurat = function(visium_slide,
                                      intra_assay, 
                                      intra_features = NULL,
                                      para_assay,
                                      para_features = NULL,
                                      l,
                                      spot_ids = NULL,
                                      out_alias = "default"){
  
  # Getting data ready to create views
  geometry = visium_slide@images$slice1@coordinates[,c(2,3)]
  
  # Intrinsic data
  intra_df = as.matrix(visium_slide@assays[[intra_assay]]@data)
  intra_df = intra_df %>%
    t %>% data.frame(check.names = F)
  
  intra_df = intra_df[rownames(geometry),]
  
  # Para data
  para_df = as.matrix(visium_slide@assays[[para_assay]]@data)
  para_df = para_df %>%
    t %>% data.frame(check.names = F)
  para_df = para_df[rownames(geometry),]
  
  
  para_pipeline(intra_df = intra_df,
                intra_features = intra_features,
                para_df = para_df,
                para_features = para_features,
                geometry = geometry,
                l = l,
                spot_ids = spot_ids,
                out_alias = out_alias)
  
}














