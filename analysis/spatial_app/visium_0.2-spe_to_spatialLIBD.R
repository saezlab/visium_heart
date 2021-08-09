# Copyright (c) [2021] [Ricardo O. Ramirez Flores, Jesus Velez]
# roramirezf@uni-heidelberg.de

#' From a list of SpatialExperiment objects, we create a spatialLIBD ready object

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)

# Load data ---------------------------------------------------------------
data_dir <- "/Users/ricardoramirez/Dropbox/PhD/Research/mi_apps/spatial_app_data/"
spe <- fs::path(data_dir, "heart_spe.rds") %>%
  readRDS()

# Add check for gene data -------------------------------------------------
rowData(spe)$gene_id <- rownames(rowData(spe))
rowData(spe)$gene_name <- rowData(spe)$symbol
rowData(spe)$gene_search <-
  paste0(rowData(spe)$gene_name, "; ", rowData(spe)$gene_id)

# Add check of QC data ----------------------------------------------------
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

# Add check for information about spot coordinates ------------------------
spatialData(spe)$barcode <- colData(spe)$Barcode

# Add check for information about samples ---------------------------------
colData(spe)$key <-
  paste0(colData(spe)$sample_id, "_", colnames(spe))

colData(spe)$Layer <- "NA"
spe$ManualAnnotation <- spe$opt_clust_integrated
colData(spe)[["clust_n"]] <- colData(spe)[["opt_clust_integrated"]] %>% strsplit(., "_") %>% map_chr(., ~ .x[2]) %>% as.numeric()

check_spe(spe)

# Save spe ----------------------------------------------------------------
data_dir %>%
  fs::path("spe_heart_spatialLIBD", ext = "rds") %>%
  saveRDS(object = spe, file = .)

