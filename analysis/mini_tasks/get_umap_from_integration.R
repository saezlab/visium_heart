# Quick script to get UMAP projections

library(tidyverse)
library(Seurat)

integrated_data <- readRDS("/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides.rds")

reductions_list <-  list(meta_data = integrated_data@meta.data,
                         reduction = integrated_data@reductions[["umap_harmony"]])

saveRDS(reductions_list,
        file = "/beegfs/work/hd_wh241/MI_revisions/processed_visium/integration/integrated_slides_umap.rds")







