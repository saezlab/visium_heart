# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Catalog of MISTy utilities
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
  
  mistyR::clear_cache()
  
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
  print(assay)
  data <- GetAssayData(visium.slide, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
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
  
  mistyR::clear_cache()
  
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      rownames_to_column() %>%
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
  run_misty(pline.views, out.alias, cached = FALSE)
}

#
# Bug in collecting results
#

collect_results_v2 <- function(folders){
  samples <- R.utils::getAbsolutePath(folders)
  message("\nCollecting improvements")
  improvements <- samples %>% furrr::future_map_dfr(function(sample) {
    performance <- readr::read_table2(paste0(sample, .Platform$file.sep, 
                                             "performance.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>% 
      dplyr::distinct()
    performance %>% dplyr::mutate(sample = sample, gain.RMSE = 100 * 
                                    (.data$intra.RMSE - .data$multi.RMSE)/.data$intra.RMSE, 
                                  gain.R2 = 100 * (.data$multi.R2 - .data$intra.R2), 
    )
  }, .progress = TRUE) %>% tidyr::pivot_longer(-c(.data$sample, 
                                                  .data$target), names_to = "measure")
  message("\nCollecting contributions")
  contributions <- samples %>% furrr::future_map_dfr(function(sample) {
    coefficients <- readr::read_table2(paste0(sample, .Platform$file.sep, 
                                              "coefficients.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>% 
      dplyr::distinct()
    coefficients %>% dplyr::mutate(sample = sample, .after = "target") %>% 
      tidyr::pivot_longer(cols = -c(.data$sample, .data$target), 
                          names_to = "view")
  }, .progress = TRUE)
  improvements.stats <- improvements %>% dplyr::filter(!stringr::str_starts(.data$measure, 
                                                                            "p\\.")) %>% dplyr::group_by(.data$target, .data$measure) %>% 
    dplyr::summarise(mean = mean(.data$value), sd = stats::sd(.data$value), 
                     cv = .data$sd/.data$mean, .groups = "drop")
  contributions.stats <- dplyr::inner_join((contributions %>% 
                                              dplyr::filter(!stringr::str_starts(.data$view, "p\\.") & 
                                                              .data$view != "intercept") %>% dplyr::group_by(.data$target, 
                                                                                                             .data$view) %>% dplyr::summarise(mean = mean(.data$value), 
                                                                                                                                              .groups = "drop_last") %>% dplyr::mutate(fraction = abs(.data$mean)/sum(abs(.data$mean))) %>% 
                                              dplyr::ungroup()), (contributions %>% dplyr::filter(stringr::str_starts(.data$view, 
                                                                                                                      "p\\.") & !stringr::str_detect(.data$view, "intercept")) %>% 
                                                                    dplyr::group_by(.data$target, .data$view) %>% dplyr::mutate(view = stringr::str_remove(.data$view, 
                                                                                                                                                           "^p\\.")) %>% dplyr::summarise(p.mean = mean(.data$value), 
                                                                                                                                                                                          p.sd = stats::sd(.data$value), .groups = "drop")), by = c("target", 
                                                                                                                                                                                                                                                    "view"))
  message("\nCollecting importances")
  importances <- samples %>% furrr::future_map(function(sample) {
    targets <- contributions.stats %>% dplyr::pull(.data$target) %>% 
      unique() %>% sort()
    views <- contributions.stats %>% dplyr::pull(.data$view) %>% 
      unique()
    maps <- views %>% furrr::future_map(function(view) {
      all.importances <- targets %>% purrr::map(~readr::read_csv(paste0(sample, 
                                                                        .Platform$file.sep, "importances_", .x, "_", 
                                                                        view, ".txt"), col_types = readr::cols()) %>% 
                                                  dplyr::distinct() %>% dplyr::rename(feature = target))
      features <- all.importances %>% purrr::map(~.x$feature) %>% 
        unlist() %>% unique() %>% sort()
      pvalues <- contributions %>% dplyr::filter(sample == 
                                                   !!sample, view == paste0("p.", !!view)) %>% dplyr::mutate(value = 1 - 
                                                                                                               .data$value)
      all.importances %>% purrr::imap_dfc(~tibble::tibble(feature = features, 
                                                          zero.imp = 0) %>% dplyr::left_join(.x, by = "feature") %>% 
                                            dplyr::arrange(.data$feature) %>% dplyr::mutate(imp = scale(.data$imp)[, 
                                                                                                                   1], `:=`(!!targets[.y], .data$zero.imp + (.data$imp * 
                                                                                                                                                               (pvalues %>% dplyr::filter(target == targets[.y]) %>% 
                                                                                                                                                                  dplyr::pull(.data$value))))) %>% dplyr::select(targets[.y])) %>% 
        dplyr::mutate(Predictor = features)
    }) %>% `names<-`(views)
  }, .progress = TRUE) %>% `names<-`(samples)
  message("\nAggregating")
  importances.aggregated <- importances %>% purrr::reduce(function(acc, 
                                                                   l) {
    acc %>% purrr::map2(l, ~(((.x %>% dplyr::select(-.data$Predictor)) + 
                                (.y %>% dplyr::select(-.data$Predictor))) %>% dplyr::mutate(Predictor = .x %>% 
                                                                                              dplyr::pull(.data$Predictor))))
  }) %>% purrr::map(~.x %>% dplyr::mutate_if(is.numeric, ~./length(samples)))
  return(list(improvements = improvements, improvements.stats = improvements.stats, 
              contributions = contributions, contributions.stats = contributions.stats, 
              importances = importances, importances.aggregated = importances.aggregated))
}











