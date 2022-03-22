# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Copying cell2location models to simpler paths
library(tidyverse)
deconv_mats_folder <- "./results/deconvolution_models/location_models/density_tables/"
deconv_rds_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"
mats_files <- list.files(deconv_mats_folder)
mats_files <- mats_files[grepl("q05", mats_files)]
rds_files <- paste0(map_chr(strsplit(mats_files,"_W"), 
                            ~ .x[1]),
                    ".rds")

walk2(mats_files, rds_files, function(mf, rf) { 
  print(rf)
  
  prefix <- paste0(gsub("[.]rds", "",rf), "_")
  
  mat <- read.csv(paste0(deconv_mats_folder, mf),
                  row.names = 1)
  
  rownames(mat) <- map_chr(strsplit(rownames(mat), prefix), ~ .x[2])
  
  colnames(mat) <- gsub("q05_spot_factors", "", colnames(mat))
  
  saveRDS(as.matrix(mat), file = paste0(deconv_rds_folder, rf))
})






