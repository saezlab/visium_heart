# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Integrate selection of slides
#' 
library(tidyverse)
library(Seurat)
library(optparse)

options(future.globals.maxSize =  4800 * 1024 ^ 2)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--all_slides"), action = "store_true", 
              default = TRUE, type = 'logical',
              help = "include all samples"),
  make_option(c("--selection"), action ="store", 
              default = "", 
              type = 'character',
              help = "samples separated by underscore"),
  make_option(c("--seed_path"), action = "store", 
              default = "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_visium/", 
              type = 'character',
              help = "where is the directory with all objects (please end with /)"),
  make_option(c("--out_path"), action= "store", 
              default = "/net/data.isilon/ag-saez/bq_shared/scellMI/processed_visium/integrated_slides.rds", 
              type = 'character',
              help = "where to save the rds object")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Parse the slides ---------------------------------------------------------------------------------
if(all_slides) {
  dsets <- c("157771", "157772",  "157775", "157777", 
             "157779", "157781", "157782", "157785")
} else {
  dsets <- unlist(strsplit(selection, "_"))
}

# Make list of data sets ----------------------------------------------------
dsets <- map(dsets, function(x) {
  
  slide_file <- paste0(seed_path,x,".rds")
  
  slide <- readRDS(slide_file)
  
  DefaultAssay(slide) <- "SCT"
  
  slide$orig.ident <- x
  
  return(slide)
  
})

# Feature selection -------------------------------------------------------------
dsets_features <- SelectIntegrationFeatures(object.list = dsets,
                                            nfeatures = 3000)

dsets_list <- PrepSCTIntegration(object.list = dsets, 
                                 anchor.features = dsets_features,
                                 verbose = FALSE)

# Integration -------------------------------------------------------------

dsets_anchors <- FindIntegrationAnchors(object.list = dsets_list, 
                                        normalization.method = "SCT", 
                                        anchor.features = dsets_features, 
                                        verbose = FALSE)

dsets_integrated <- IntegrateData(anchorset = dsets_anchors, 
                                  normalization.method = "SCT", 
                                  verbose = FALSE)

# Analysis of integration -------------------------------------------------------------

DefaultAssay(dsets_integrated) <- "integrated"
dsets_integrated <- RunPCA(dsets_integrated, verbose = FALSE)
dsets_integrated <- RunUMAP(dsets_integrated, dims = 1:30)
dsets_integrated <- FindNeighbors(dsets_integrated, 
                                  reduction = "pca", 
                                  dims = 1:30,
                                  k.param = 5)
dsets_integrated <- FindClusters(dsets_integrated, 
                                 resolution = 1)

# Saving objects -------------------------------------------------------------

saveRDS(dsets_integrated, file = out_path)

meta_data <- FetchData(dsets_integrated,
                       vars = c("orig.ident",
                                "lt_id",
                                "UMAP_1", 
                                "UMAP_2")) %>%
  as.data.frame()

out_path_meta <- gsub(pattern = "[.]rds",
                      replacement = "_metadata.rds",
                      out_path)

saveRDS(meta_data, file = out_path_meta)
