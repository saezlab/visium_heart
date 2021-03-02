# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Estimate colocalization with ct scores in
#' all slides
#' 
#' Fits NMF in a unique matrix
#' 
#' 
library(optparse)
library(Seurat)
library(tidyverse)
library(NNLM)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--deconv_mats_folder"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "location where deconvolution matrices are located (rds), features in columns"),
  make_option(c("--slide_files_folder"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where do the visium rds objects are located?"),
  make_option(c("--out_slides_folder"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the analysis objects"),
  make_option(c("--alias"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "assay alias"),
  make_option(c("--ks"), 
              action= "store", 
              default = "6,8,10,12", 
              type = 'character',
              help = "number of factors in NMF, separated by commas")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
slide_files <- list.files(slide_files_folder)

data_df <- tibble(slide_name = gsub("[.rds]","",slide_files),
                  ct_mat = paste0(deconv_mats_folder, sort(deconv_mats)),
                  slide_file = paste0(slide_files_folder, sort(slide_files))) %>%
  dplyr::mutate(slides = map(slide_file, readRDS),
                deconv_mats = map(ct_mat, readRDS)) %>%
  dplyr::select(slide_name, slides, deconv_mats) %>%
  dplyr::mutate(slides = map2(slides, slide_name, function(x, y) {
    x$sample <- y
    return(x)
  })) %>%
  dplyr::mutate(slides = map2(slides, deconv_mats, function(x, y) {
    x[[alias]] <- CreateAssayObject(data = t(y))
    return(x)
  }))

integrated_slide <- reduce(data_df$slides,
                 merge,
                 merge.data = TRUE)

rm(data_df)

# Now run the NMF --------------------------------------

ks <- as.numeric(unlist(strsplit(x = ks, split = ",")))

walk(ks, function(k) {
  
  print(k)
  
  decomp <- nnmf(t(as.matrix(integrated_slide@assays[[alias]]@data)), 
                 k = k, 
                 max.iter = 500)
  
  # Making objects compatible with Seurat structure
  factor_weights <- decomp$W
  colnames(factor_weights) <- paste0("factor_",1:ncol(factor_weights))
  
  factor_loadings <- decomp$H
  rownames(factor_loadings) <- paste0("factor_",1:nrow(factor_loadings))
  
  NMF_res <- list("factor_weights" = factor_weights,
                  "factor_loadings" = factor_loadings,
                  "meta_data" = integrated_slide@meta.data)
  
  saveRDS(file = paste0(out_slides_folder,"_k",k,".rds"), 
          NMF_res)
  
})

