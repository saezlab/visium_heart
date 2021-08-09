# Copyright (c) [2021] [Ricardo O. Ramirez Flores, Jesus Velez]
# roramirezf@uni-heidelberg.de

#' Converts a collection of visium seurat datasets into a list of SpatialExperiments

# Load libraries ----------------------------------------------------------
library(dplyr)
library(purrr)
library(SpatialExperiment)
library(Seurat)
library(SeuratObject)
library(optparse)


# Command line parameters -------------------------------------------------
option_list <- list(
  make_option(c("--raw_folder"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%sample(/outs)/filtered_feature_bc_matrix"),
  make_option(c("--processed_folder"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "is the path added a folder with structure ./%sample(/outs)/filtered_feature_bc_matrix"),
  make_option(c("--out_folder"), 
              action ="store", 
              default = "default", 
              type = 'character',
              help = "output folder"),
  make_option(c("--outs_structure"), 
              action = "store", 
              default = TRUE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%sample/outs/filtered_feature_bc_matrix")
  )

# Specify working directories ---------------------------------------------
visium_data_dir <- raw_folder
seurat_data_dir <- processed_folder
seurat_data_dir <- fs::path(data_dir, "objects")

# Small function to add outs structure if necessary -----------------------
paste_outs <- function(path, outs_structure) {
  if(outs_structure) {
    paste0(path, "/outs")
  } else {
    path
  }
  
}

# fs::dir_tree(data_dir, type = "dir", recurse = FALSE)
# data/raw
# ├── processed_visium_revisions
# └── visium_data

# Load visiums as Spatial Experiments -------------------------------------
visium_spe <- visium_data_dir %>% 
  fs::dir_ls(type = "dir") %>% 
  purrr::set_names(nm = fs::path_file(.)) %>% 
  purrr::map_chr(., paste_outs, outs_structure) %>%
  SpatialExperiment::read10xVisium() %>% 
  base::identity()

# Load seurat objects and merge them --------------------------------------
merged_seurat2 <- seurat_data_dir %>%
  fs::dir_ls(glob = "*.rds") %>% 
  rlang::set_names(fs::path_ext_remove(fs::path_file(.))) %>% 
  purrr::map(~{
    seurat_object <- readRDS(.x)
    seurat_object
  }) %>% 
  {
    # Seurat:::merge.SCTAssay(
    merge(
      x = .[[1]],
      y = .[2:length(.)],
      add.cell.ids = names(.)
    )
  } %>%
  base::identity()

# Copy seurat metadata into visium_spe ------------------------------------
merged_seurat <- merged_seurat2
barcodes_ok <- visium_spe %>% 
  {
    all(
      paste0(colData(.)$sample_id, "_", rownames(colData(.))) ==
        colnames(merged_seurat)
    )
  }    

get_seurat_assays_data <- function(seurat) {
  seurat@assays %>% 
    purrr::list_modify(
      Spatial = purrr::zap(),
      SCT = purrr::zap()
    ) %>% 
    purrr::discard(purrr::is_empty) %>%
    purrr::map2(
      .y = names(.),
      .f = ~{
        .x@data %>% 
          `rownames<-`(paste0(.y, "_", rownames(.)))
      }
    ) %>%
    purrr::reduce(base::rbind) %>% 
    `colnames<-`(
      stringr::str_remove(
        string = colnames(.),
        pattern = ".+_"
      )
    ) %>% 
    t() %>% 
    base::identity()
}


if (barcodes_ok) {
  visium_seurat_spe <- visium_spe
  colData(visium_seurat_spe) <- cbind(
    colData(visium_spe),
    merged_seurat@meta.data,
    get_seurat_assays_data(merged_seurat)
  )
  
  # Save object -------------------------------------------------------------
  out_dir <- here::here("spatial_app_data") %>% 
    fs::dir_create()
  
  visium_seurat_spe %>% 
    saveRDS(file = fs::path(out_dir, "heart_spe", ext = "rds"))
}
