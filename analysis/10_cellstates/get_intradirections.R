# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate ct abundances with cell-states
#' Correlate cell-states with PROGENy
#' Correlate PROGENY with ct-abundances
#' 

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)

# Import path pointers

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()


# Get annotations of where are the actual cell-types -----------------------------------------------

cell_flags <- map(visium_df$visium_file, function(visium_file) {
  print(visium_file)
  slide <- readRDS(visium_file)
  meta <- slide@meta.data %>%
    rownames_to_column("spot_id")
  
  flag_ids <- grepl("flag",colnames(meta))
  flag_ids <- colnames(meta)[flag_ids]
  
  meta[,c("spot_id", flag_ids)]
})

names(cell_flags) <- visium_samples
  
cell_flags <- enframe(cell_flags,name = "sample_id") %>%
  unnest() %>%
  mutate(general_id = paste0(sample_id, "..",spot_id))


# Get a joint matrix with cell scores ----------------------------------------------------------

get_assay_matrices <- function(visium_file_vect, visium_samples_vect, assay) {
  # Get mat for each file
  mat_list <- map2(visium_file_vect, visium_samples_vect,function(visium_file, visium_sample) {
    print(visium_sample)
    slide <- readRDS(visium_file)
    mat <- GetAssayData(slide, assay = assay) %>% as.matrix() %>% t()
    rownames(mat) <- paste0(visium_sample, "..", rownames(mat))
    return(mat)
  })
  
  # Put them together
  names(mat_list) <- paste0("sample", visium_samples_vect)
  mat <- purrr::reduce(mat_list, rbind)
}

# c2l scores  -----------------------------------------------
cell_props <- get_assay_matrices(visium_file_vect = visium_df$visium_file,
                                 visium_samples_vect = visium_df$sample,
                                 assay = "c2l")
  
# Get a joint matrix with progeny scores -----------------------------------------------
cell_paths <- get_assay_matrices(visium_file_vect = visium_df$visium_file,
                                 visium_samples_vect = visium_df$sample,
                                 assay = "progeny")

# Get a joint matrix with cell states -----------------------------------------------
cell_states <- get_assay_matrices(visium_file_vect = visium_df$visium_file,
                                  visium_samples_vect = visium_df$sample,
                                  assay = "cell_states")

# Create same cell order:
cell_order <- cell_flags$general_id
cell_props <- cell_props[cell_order, ]
cell_paths <- cell_paths[cell_order, ]
cell_states <- cell_states[cell_order, ]

# For each cell we will only select spots of interest: flag = 1

get_ct_corrs <- function(cell_flags, ct) {
  
  useful_ids <- cell_flags[[paste0(ct, "_flag")]] == 1
  useful_ids <- cell_flags[["general_id"]][useful_ids]
  
  red_cell_props <- cell_props[, !grepl(ct, colnames(cell_props))]
  #red_cell_props_A <- cell_props
  red_states <- cell_states[, grepl(ct, colnames(cell_states))]
  all_states <- cell_states
  
  list(ct_ct = cor(red_cell_props[useful_ids, ], red_cell_props[useful_ids, ]),
       ct_paths = cor(red_cell_props[useful_ids, ], cell_paths[useful_ids, ]),
       ct_states = cor(red_cell_props[useful_ids, ], red_states[useful_ids, ]),
       states_paths = cor(cell_paths[useful_ids, ], red_states[useful_ids, ]),
       states = cor(red_states[useful_ids, ], red_states[useful_ids, ]),
       states_all = cor(all_states[useful_ids, ], all_states[useful_ids, ]))
  
}

walk(c("Fib", "Myeloid"), function(ct) {
  
  print(ct)
  
  list_res <- get_ct_corrs(cell_flags, ct = ct)
  
  pdf_out <- paste0("./results/state_structure/",ct,"/correlation_summary.pdf")
  
  pdf(pdf_out, height = 4, width = 5)
  
    walk(list_res, function(x) {
    draw(ComplexHeatmap::Heatmap(x,))
    })
  
  dev.off()
  
})

#cardio_res <- get_ct_corrs(cell_flags, ct = "cardiomyocyte")
fibro_res <- get_ct_corrs(cell_flags, ct = "Fib")
#endo_res <- get_ct_corrs(cell_flags, ct = "endothelial")

