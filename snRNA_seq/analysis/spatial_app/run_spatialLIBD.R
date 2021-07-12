# Copyright (c) [2021] [Ricardo O. Ramirez Flores, Jesus Velez]
# roramirezf@uni-heidelberg.de

#' Deploy spatialLIBD app

library(SpatialExperiment)
library(spatialLIBD)

# Load data ---------------------------------------------------------------
data_dir <- "/Users/ricardoramirez/Dropbox/PhD/Research/mi_apps/spatial_app_data/"
spe_heart <- data_dir %>%
  fs::path("spe_heart_spatialLIBD", ext = "rds") %>%
  readRDS()

variables <- colnames(colData(spe_heart))
numeric_cols <- sapply(colData(spe_heart), is.numeric)

spe_continuous_vars <- variables[numeric_cols]
spe_discrete_vars <- variables[!numeric_cols]

spe_discrete_vars <- setdiff(
  spe_discrete_vars,
  c(
    "Sample",
    "Barcode",
    "sample_id",
    "key",
    "Layer",
    "orig.ident",
    "ManualAnnotation"
  )
)

colData(spe_heart)[["opt_clust_integrated"]] <- colData(spe)[["opt_clust_integrated"]] %>% 
  strsplit(., "_") %>% 
  map_chr(., ~ .x[2]) %>% 
  as.numeric(.)

colData(spe_heart)[["opt_clust_integrated"]] <- colData(spe_heart)[["opt_clust_integrated"]]  + 1


run_app(
  spe = spe_heart,
  sce_layer = NULL,
  modeling_results = NULL,
  sig_genes = NULL,
  spe_discrete_vars = spe_discrete_vars,
  spe_continuous_vars = spe_continuous_vars,
  default_cluster = "opt_clust_integrated"
)









