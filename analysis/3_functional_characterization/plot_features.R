# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we explore markers from an integrated 
#' data set

library(Seurat)
library(tidyverse)
library(cowplot)
library(optparse)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--assay_name"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--features"), 
              action ="store", 
              default = "all", 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--slot"), 
              action ="store", 
              default = "data", 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--group_by"), 
              action ="store", 
              default = "data", 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--out_fig_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the fig objects"))

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Main ------------------------------------------------------------------

visium_data <- readRDS(path)

if(features == "all") {
  
  features <- rownames(GetAssayData(visium_data,
                                    slot = slot, 
                                    assay = assay_name))
  
} else {
  
  features <- unlist(strsplit(features, ","))
  
}

violin_list <- VlnPlot(visium_data, features = features, 
                       group.by = group_by,
                       assay = assay_name, 
                       slot = slot, combine = F)

DefaultAssay(visium_data) <- assay_name

umap_ps <- FeaturePlot(visium_data, 
                       features = features,
                       combine = F, raster = T)

pdf(out_fig_path, height = 7, width = 7)

walk(violin_list, print)
walk(umap_ps, print)

dev.off()

