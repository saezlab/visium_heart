# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we will run MISTy to explain the cell-type scores of different cell-types
#' First, we require to define meta variables per spot that flag spots with
#' enough cell-type
#' 
#' Later, we used those unified tags to get specific data as shown in niche_misty.R
#' and run MISTy pipelines

library(Seurat)
library(mistyR)
library(tidyverse)
source("./analysis/utils/misty_pplne_utils.R")


# Import path pointers

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()


# First for each slide we will create metavariables that flag the location of a cell-type in a spot

c2l_assay <- "c2l_major_props"

# Eveything except cardiomyocytes and fibroblast must represent 10% of celltype score
ct_prop_param <- tibble(cts = c("adipocytes", "cardiomyocytes", "endothelial-cells", "fibroblasts",
                                "lymphatic-endo", "macrophages", "neuronal", "pericytes", "t-cells")) %>%
  mutate(prop_param = ifelse(cts %in% c("cardiomyocytes", "fibroblasts"), 0.1, 0.1))

add_ct_flags <- function(slide, ct_prop_param, cell_props) {
  
  for(ct in ct_prop_param$cts) {
    
    ix <- grepl(ct, ct_prop_param$cts)
    
    prop_param <- ct_prop_param$prop_param[ix]
    
    slide[[paste0(ct,"_flag")]] = ifelse(cell_props[, ct] <= prop_param, 0, 1)
    
  }
  
  return(slide)
  
}

walk(visium_df$visium_file, function(visium_file) {
  
  print(visium_file)
  
  slide <- readRDS(visium_file)
  cell_props <- GetAssayData(slide, assay = c2l_assay) %>% t()
  
  slide_ct_prop_param <- ct_prop_param %>%
    dplyr::filter(cts %in% colnames(cell_props))
  
  slide <- add_ct_flags(slide = slide,
                        ct_prop_param = slide_ct_prop_param,
                        cell_props = cell_props)
  
  saveRDS(slide, file = visium_file)
  
})

# Once each slide has been modified to have cell-type flags we fit a niche MISTy analysis

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
                                    spot_group,
                                    group_id,
                                    view_types = list("main" = "intra", 
                                                      "juxta" = "juxta",
                                                      "para" = "para"),
                                    view_params = list("main" = NULL, 
                                                       "juxta" = 5,
                                                       "para" = 15),
                                    view_assays = list("main" = "SCT",
                                                        "juxta" = "SCT",
                                                        "para" = "SCT"),
                                    feature_exception = list("main" = NULL,
                                                             "juxta" = NULL,
                                                             "para" = NULL)) {
  
  print(slide_path)
  # Reading data
  visium_slide <- readRDS(slide_path)
  
  # This identifies the specific spots to be used
  Idents(visium_slide) <- spot_group
  spots <- WhichCells(visium_slide, idents = group_id)
  
  # Get useful features by default all the features in the assay are included
  view_features <- map2(view_assays, feature_exception, function(assay_name, exception) {
    
    cell_feats <- rownames(visium_slide@assays[[assay_name]])
    
    if(!is.null(exception)) {
      cell_feats <- cell_feats[!cell_feats %in% exception]
    }
    
    return(cell_feats)
    
  })
  
  get_niche_views(visium.slide = visium_slide,
                  view.assays = view_assays,
                  view.features = view_features,
                  view.types = view_types,
                  view.params = view_params,
                  spot.ids = spots)
  
}

# Main pipeline definition

# Get niche relevant data

# In this first step we subset data by a variable of interest

view_names <- c("main_state", "intra_props", 
                "juxta_props", "para_props") %>%
  set_names()

state_vector <- c("cardiomyocytes-0", "cardiomyocytes-1", "cardiomyocytes-2", "cardiomyocytes-3",
                  "endothelial-cells-0", "endothelial-cells-1", "endothelial-cells-2", "endothelial-cells-3", "endothelial-cells-4", "endothelial-cells-5",
                  "fibroblasts-0", "fibroblasts-1", "fibroblasts-2", "fibroblasts-3", "fibroblasts-4", "fibroblasts-5",
                  "macrophages-0", "macrophages-1",  "macrophages-2", "macrophages-3",
                  "macrophages-4", "macrophages-5", "macrophages-6", "macrophages-7",
                  "pericytes-0", "pericytes-1", "pericytes-2", "pericytes-3")

state_misty_pplne <- function(state_flag, 
                              visium_df,
                              state_vector,
                              out_alias = "./visium_results_manuscript/cell_state_misty/") {
  
  print(state_flag)
  
  major_ct <- state_flag %>% gsub("_flag", "", .) %>% gsub("[.]","-", .)
  ct_exception_ix <- grepl(pattern = major_ct, x = state_vector)
  ct_exception <- state_vector[!ct_exception_ix]
  
  # Get view data from each slide for the selected cell-type flag
  # filtering states to only be daughters of the cell-type
  slide_views <- visium_df %>%
    mutate(niche_data = map(visium_file, 
                            get_blocked_misty_views,
                            spot_group = state_flag,
                            group_id = 1,
                            # Define spatial context of each view -----
                            view_types = list("main_state" = "intra", 
                                              "intra_props" = "intra",
                                              "juxta_props" = "juxta",
                                              "para_props" = "para"),
                            # Define additional parameters (l in case of paraview,
                            # n of neighbors in case of juxta) --------
                            view_params = list("main_state" = NULL, 
                                               "intra_props" = NULL,
                                               "juxta_props" = 5,
                                               "para_props" = 15),
                            view_assays = list("main_state" = "cell_states",
                                               "intra_props" = "c2l_major_props",
                                               "juxta_props" = "c2l_major_props",
                                               "para_props" = "c2l_major_props"),
                            feature_exception = list("main_state" = ct_exception,
                                                     "intra_props" = c(major_ct, "adipocytes"),
                                                     "juxta_props" = c(major_ct, "adipocytes"),
                                                     "para_props" = c(major_ct, "adipocytes"))))
  
  # Second step is to create alternative joint views:
  # First per view name you extract data
  
  niche_dat <- slide_views %>%
    mutate(niche_data = map(niche_data, function(slide_dat) {
      
      view_dat <- map(view_names, function(vn) { 
        
        print(vn)
        first_level <- slide_dat[[vn]]
        second_level_names <- names(first_level)
        second_level_ix <- second_level_names[grepl(vn, second_level_names)]
        second_level <- first_level[[second_level_ix]]
        second_level[["data"]]
          
        })
    })) %>%
    unnest_wider("niche_data")
  
  # Third step is to format the joint views for MISTy
  
  niche_dat <- map(view_names, ~ bind_rows(niche_dat[[.x]]))
  
  views_main <- create_initial_view(niche_dat[[view_names[1]]], unique.id = NULL)
  
  other_views <- map(view_names[-1], function(x) {
    
    create_view(name = x, 
                data = niche_dat[[x]],
                abbrev = x)
    
  })
  
  for(view in view_names[-1]) {
    print(view)
    views_main <- add_views(
      views_main,
      other_views[[view]]
    )
    
  }
  
  # Define path
  misty_out <- paste0(out_alias, major_ct)
  
  run_misty(views_main, results.folder = misty_out, cached = FALSE)
  
  
}


cell_flags <- set_names(c("cardiomyocytes_flag", "fibroblasts_flag", 
                          "pericytes_flag", "macrophages_flag", 
                          "endothelial.cells_flag"))

test <- map(cell_flags, 
            state_misty_pplne,  
            visium_df = visium_df, 
            state_vector = state_vector,
            out_alias = "./visium_results_manuscript/cell_state_misty/props_")


walk(test, function(misty_out_path) {
  
  pdf_file <- paste0(misty_out_path, ".pdf")
  
  misty_out <- mistyR::collect_results(misty_out_path)
  
  pdf(file = pdf_file, height = 5, width = 5)
  
  print(misty_out$improvements %>%
          dplyr::filter(grepl("R2", measure) &
                          !grepl("p", measure) &
                          !grepl("gain", measure)) %>%
          ggplot(aes(x = target, fill = measure, y = value)) + 
          geom_bar(stat="identity",position = "dodge") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  
  misty_out %>% mistyR::plot_improvement_stats()
  misty_out %>% mistyR::plot_view_contributions()
  
  misty_out %>% mistyR::plot_interaction_heatmap("intra", cutoff = 0.5) 
  
  misty_out %>% mistyR::plot_interaction_heatmap("intra_props", cutoff = 0.5)
  
  misty_out %>% mistyR::plot_interaction_heatmap("juxta_props", cutoff = 0.5)
  
  misty_out %>% mistyR::plot_interaction_heatmap("para_props",cutoff = 0.5)
  
  dev.off()
  
})



























