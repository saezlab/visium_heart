# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we will generate subsamples of the single cell atlas
#' This is to run cell2location models with a better estimate

library(SingleCellExperiment)
library(tidyverse)
library(HDF5Array)
library(zellkonverter)

# Set magic seed
set.seed(241099)
sampling_perc <- 0.4
sc_data <- loadHDF5SummarizedExperiment("./processed_snrnaseq/integration/integrated_rnasamples_ann_sce/")
out_file <- "./processed_snrnaseq/subsamples/ssample_integrated_rnasamples_ann"

# Get number of cells per cell-type
group_var <- "cell_type"

cell_pool <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("cell_ID") %>%
  group_by_at(group_var) %>%
  dplyr::select(cell_ID) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  dplyr::rename("cell_ids" = data) %>%
  mutate(ncells = map_dbl(cell_ids, length)) %>%
  mutate(samplecells = floor(ncells * sampling_perc))

subsample_atlas <- function(sc_data, cell_pool) {
  
  subsample_cells <- map2(cell_pool$cell_ids, 
                          cell_pool$samplecells, sample) %>%
    unlist()
  
  subsampled_atlas <- sc_data[ , subsample_cells]
  
}

seq(1:5) %>% walk(., function(iter) {
  print(iter)
  
  ss_atlas <- subsample_atlas(sc_data, cell_pool)
  
  final_out_file <- paste0(out_file, "_", iter,".h5ad")
  
  print(paste0("saving in ", final_out_file))
  
  writeH5AD(ss_atlas, file = final_out_file)
  
})





