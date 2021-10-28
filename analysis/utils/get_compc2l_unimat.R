# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script I generate a unified matrix of cell2location scores and compositions
#' I also provide utilities to create flag matrices to filter spots in a faster fashion

library(tidyverse)

# Here we will get the matrix of compositions using c2l results ----------------------------------

c2l_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"
assay_name <- "c2l"

# Get atlas meta:
atlas_meta <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["annotations"]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id,"_"), ~.x[[1]])) %>%
  dplyr::mutate(row_id = paste0(orig.ident, "..", cell_id))

# Get cell2location files --------------------------------
c2l_files <- list.files(c2l_folder, full.names = F)

c2l_samples <- map_chr(strsplit(c2l_files,".rds"), 
                       ~ .x[1])

c2l_df <- tibble(c2l_file = paste0(c2l_folder, 
                                   c2l_files),
                 sample = c2l_samples) 

# Generates list of matrix of c2l converted proportions
list_matrices <- map2(c2l_df$c2l_file, c2l_df$sample, function(f, s) {
  mat <- readRDS(f)
  rownames(mat) <- paste0(s, "..", rownames(mat))
  prop_mat <- base::apply(mat, 1, function(x) {
    
    x/sum(x)
    
  })
  
  return(t(prop_mat))
})

# Reduces this to a single matrix
names(list_matrices) <- c2l_df$sample
integrated_compositions <- purrr::reduce(list_matrices, rbind)
rm(list_matrices)
# Order it as meta data
integrated_compositions <- integrated_compositions[atlas_meta$row_id, ]

integrated_compositions <- integrated_compositions %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::mutate("sample_id" = strsplit(id, "[..]") %>%
                  map_chr(., ~ .x[[1]]),
                "spot_id" = strsplit(id, "[..]") %>%
                  map_chr(., ~ .x[[3]])) %>%
  dplyr::select(-id) %>%
  pivot_longer(-c("sample_id", "spot_id"))

saveRDS(integrated_compositions, file = "./markers/comp2cl.rds")
