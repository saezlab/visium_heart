# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Map shared factors - In theory you already have the counts in the slide objects

all_nmf_folder <- "./visium_results_manuscript/colocalization/"
nmf_res <- list.files(all_nmf_folder)[grepl("[.]rds",list.files(all_nmf_folder))]
nmf_id <- set_names(gsub("[.]rds", "", nmf_res))
  
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)

assign_sharedfactors <- function(slide_file, nmf_res_file, alias_nmf) {
  print(slide_file)
  # Read spatial transcriptomics data
  slide_id <- gsub("[.]rds", "", slide_file)
  slide <- readRDS(paste0(slide_files_folder, slide_file))
  spot_order <- map_chr(strsplit(colnames(slide),"-"), ~.x[1])
  
  # Read NMF results ------------------------------------------------
  nmf_res_list <- readRDS(paste0(all_nmf_folder, nmf_res_file))
  
  # Reading meta_data and filtering for useful spots, both factors and meta
  
  meta_data <- nmf_res_list$meta_data %>%
    rownames_to_column("cell_id_raw") %>%
    dplyr::mutate("cell_id" = map_chr(strsplit(cell_id_raw,"-"), ~.x[1])) %>%
    dplyr::filter(sample == slide_id) %>%
    dplyr::select(sample,cell_id_raw, cell_id)
  
  factor_weights <- nmf_res_list$factor_weights[meta_data$cell_id_raw,]
    
  # Order data to be in line with original slide
  rownames(factor_weights) <- meta_data$cell_id
  factor_weights <- factor_weights[spot_order,]
  rownames(factor_weights) <- colnames(slide)
  
  # Add assay
  slide[[alias_nmf]] <- CreateAssayObject(data = t(factor_weights))
  
  DefaultAssay(slide) <- alias_nmf
  
  all_factors_plts <- SpatialFeaturePlot(slide, features = rownames(slide), ncol = 3)
  
  return(all_factors_plts)
}

k_res <- walk(nmf_id, function(x){ 
  
  print(x)
  nmf_res_file <- paste0(x,".rds")
  
  factor_map_plts <- map(slide_files, 
      assign_sharedfactors, 
      alias_nmf = x, 
      nmf_res_file = nmf_res_file)
  
  # Read NMF results ------------------------------------------------
  nmf_res_list <- readRDS(paste0(all_nmf_folder, nmf_res_file))
  
  # Plotting loadings ------------------------------------------------
  factor_loadings <- nmf_res_list$factor_loadings %>%
    as.data.frame() %>%
    rownames_to_column("factor") %>%
    pivot_longer(-factor, 
                 names_to = "feature", 
                 values_to = "loading")
  
  loadingplts <- ggplot(factor_loadings, aes(fill = loading, x = factor, y = feature)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  pdf(file = paste0(all_nmf_folder,x,".pdf"), height = 15, width = 15)
  
  print(loadingplts)
  
  walk(factor_map_plts, print)
  
  dev.off()
  
})















