# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we prototype niche comparison
#' First we select certain spots per slide by a given feature (niche definition)
#' Then we create frankenstein views that are fitted to a single MISTy model
#' 
#' 
#' 

library(tidyverse)
library(Seurat)
library(scater)
library(mistyR)
source("./analysis/utils/misty_pplne_utils.R")

# data collection functions
get_niche_views <- function(visium.slide,
                            view.assays,
                            view.features = NULL,
                            view.types,
                            view.params,
                            spot.ids = NULL) {
  
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  get_misty_views(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids
  )
  
}

get_blocked_misty_views <- function(slide_path,
                                    assay_name,
                                    spot_group,
                                    group_id,
                                    # Define spatial context of each view -----
                                    view_types = list("main" = "intra", 
                                                       "juxta" = "juxta",
                                                       "para" = "para"),
                                    # Define additional parameters (l in case of paraview,
                                    # n of neighbors in case of juxta) --------
                                    view_params = list("main" = NULL, 
                                                        "juxta" = 5,
                                                        "para" = 15),
                                    feature_exception = NULL) {
  
  print(slide_path)
  # Reading data
  visium_slide <- readRDS(slide_path)
  
  # This identifies the specific spots to be used
  Idents(visium_slide) <- spot_group
  spots <- WhichCells(visium_slide, idents = group_id)
  
  # Get useful features by default all the features in the assay are included
  cell_feats <- rownames(visium_slide@assays[[assay_name]])
  if(!is.null(feature_exception)) {
    cell_feats <- cell_feats[!cell_feats %in% feature_exception]
  }
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay_name,
                      "juxta" = assay_name,
                      "para" = assay_name)
  
  # Define features of each view ------------
  view_features <- list("main" = cell_feats, 
                            "juxta" = cell_feats,
                            "para" = cell_feats)
      
  get_niche_views(visium.slide = visium_slide,
                      view.assays = view_assays,
                      view.features = view_features,
                      view.types = view_types,
                      view.params = view_params,
                      spot.ids = spots)
      
  }
  
# Main

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()

# Get relevant niches --------------------------------------------------

niche_gex_obj <- readRDS("./processed_visium/integration/ps_integrated_slides_niches.rds")[[1]]

meta_data <- colData(niche_gex_obj$gex) %>%
  as.data.frame() %>%
  dplyr::rename(sample = orig.ident,
                niche = opt_clust_integrated) %>%
  dplyr::mutate(niche = paste0("niche_", niche)) %>%
  dplyr::mutate(col_id = paste0(sample, ".", niche)) %>%
  dplyr::select(sample, niche)

visium_df <- visium_df %>%
  left_join(meta_data)

# Get niche relevant data

# In this first step we subset data by a variable of interest

view_names <- c("main", "juxta", "para") %>%
  set_names()

# Here you do a map over all niches of interest ----------- TO DO

niche_misty_pplne <- function(niche_selection, 
                              view_names, 
                              out_alias = "./visium_results_manuscript/niche_misty/") {
  
  print(niche_selection)
  
  # Get view data from each slide for the selected niche
  niche_dat <- visium_df %>%
    group_by(niche) %>%
    nest() %>%
    dplyr::filter(niche == niche_selection) %>%
    mutate(data = map2(niche, data, function(niche_id, df) {
      
      slide_views <- df %>%
        mutate(niche_data = map(visium_file, get_blocked_misty_views,
                                assay_name = "c2l_major_props",
                                spot_group = "opt_clust_integrated",
                                group_id = niche_id,
                                # Define spatial context of each view -----
                                view_types = list("main" = "intra", 
                                                  "juxta" = "juxta",
                                                  "para" = "para"),
                                # Define additional parameters (l in case of paraview,
                                # n of neighbors in case of juxta) --------
                                view_params = list("main" = NULL, 
                                                   "juxta" = 5,
                                                   "para" = 15),
                                feature_exception = NULL))
      
    }
    ))
  
  # Second step is to create alternative joint views:
  # First per view name you extract data
  
  niche_dat <- niche_dat %>%
    mutate(data = map(data, function(slide_dat) {
      
      view_dat <- map(view_names, function(vn) { 
        
        print(vn)
        
        map(slide_dat$niche_data, function(misty_views) {
          
          first_level <- misty_views[[vn]]
          second_level_names <- names(first_level)
          second_level_ix <- second_level_names[grepl(vn, second_level_names)]
          second_level <- first_level[[second_level_ix]]
          second_level[["data"]]
          
        })
        
      })
      
    })) %>% 
    unnest_wider("data")
  
  # Third step is to format the joint views for MISTy
  main = bind_rows(niche_dat$main)
  juxta = bind_rows(niche_dat$juxta)
  para = bind_rows(niche_dat$para)
  
  views_main <- create_initial_view(main, unique.id = NULL)
  
  other_views <- list("juxta" = create_view(name = "juxta", 
                                            data = juxta,
                                            abbrev = "juxta"),
                      "para" = create_view(name = "para", 
                                           data = para, 
                                           abbrev = "para"))
  
  views_main <- add_views(
    views_main,
    unlist(other_views, recursive = FALSE)
  )
  
  # Define path
  misty_out <- paste0(out_alias, niche_selection)
  
  run_misty(views_main, results.folder = misty_out)
  
  
}

# Main: Run 

niche_list <- c("niche_0", "niche_1", "niche_5", "niche_8") %>%
  set_names()

test <- map(niche_list, niche_misty_pplne,  
            view_names = view_names, 
            out_alias = "./visium_results_manuscript/niche_misty/")


walk(test, function(misty_out_path) {
  
  pdf_file <- paste0(misty_out_path, ".pdf")
  
  misty_out <- mistyR::collect_results(misty_out_path)
  
  pdf(file = pdf_file, height = 8, width = 8)
  
  misty_out %>% mistyR::plot_improvement_stats()
  misty_out %>% mistyR::plot_view_contributions()
  
  misty_out %>% mistyR::plot_interaction_communities("intra")
  misty_out %>% mistyR::plot_interaction_heatmap("intra", cutoff = 0.5) 
  
  misty_out %>% mistyR::plot_interaction_communities("para")
  misty_out %>% mistyR::plot_interaction_heatmap("para", cutoff = 0.5)
  
  misty_out %>% mistyR::plot_interaction_communities("juxta")
  misty_out %>% mistyR::plot_interaction_heatmap("juxta",cutoff = 0.5)
  
  dev.off()

})