# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Calculate cell-cycle phase 
library(Seurat)
library(tidyverse)

sn_dat <- readRDS("/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/integrated_rnasamples_ann.rds")
Idents(sn_dat) <- "cell_type"

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "./markers/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sn_dat <- CellCycleScoring(sn_dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cell_cycle_data <- sn_dat@meta.data %>% 
  as.data.frame() %>%
  rownames_to_column("nuclei_id") %>%
  dplyr::select_at(c("nuclei_id", "S.Score", "G2M.Score", "Phase"))

umap_data <- sn_dat@reductions$umap_harmony@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("nuclei_id")

final_obj <- list("cell_cycle" = cell_cycle_data, 
                  "umap" = umap_data)

saveRDS(final_obj, file = "/beegfs/work/hd_wh241/MI_revisions/processed_snrnaseq/integration/cell_cycle_obj.rds")
