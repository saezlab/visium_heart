# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Catalog of MISTy utilities for spatial interaction analysis
#' 

run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {
  
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
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}

# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  data <- GetAssayData(visium.slide, assay = assay) %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% dplyr::slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    
    view.data.tmp <- view.data.init
    data.ix <- "intraview"

    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
    
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
    
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  clear_cache(view.data.init$misty.uniqueid)
  
  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias)
}

# Gets view data without running pipelines
get_misty_views <- function(view.data,
                            view.features,
                            view.types,
                            view.params,
                            geometry,
                            spot.ids = NULL) {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt,
    view.types,
    view.params,
    view.names),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  # Run MISTy
  return(all.views)
}