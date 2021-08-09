# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we verify that removed spots have less nuclei counts
#' than high quality spots
#' Visium data structure
#' folder
#' |
#' sample---outs
#'          |
#'          ---spatial
#'          ---filtered_feature_bc_matrix.h5
#'          
#' 1) Read 10x files unfiltered
#' 2) Read nuclei quantification
#' 3) Correlate ngenes, ncounts, and nuclei
#' 4) Make violin plots of filtered vs unfiltered

library(tidyverse)
library(Seurat)
library(cowplot)

folder <- TRUE
path <- "./visium_data/"
nuclei_quant_path <- "./nuclei_quant/tissue_spot_counts_"
out_qc_pdf <- "./processed_visium/initial_qc/tissue_spot_counts_qc_"

# This is an option to give a folder and parse all the files to process --------------------------------
if(folder) {
  sample_names <- list.files(path)
  slide_files <- paste0(path,
                        sample_names,
                        "/outs")
} else {
  sample_names <- id
  slide_files <- path
}

# Make data frame of parameters ------------------------------------------------------------------------
param_df <- tibble(sample_name = sample_names,
                   slide_file = slide_files,
                   nuclei_file = paste0(nuclei_quant_path, sample_names, ".csv"),
                   pdf_file = paste0(out_qc_pdf, sample_names, ".pdf")) %>%
  dplyr::filter(nuclei_file %in% paste0("./nuclei_quant/", 
                                        list.files("./nuclei_quant/")))

# Automatic processing of each data set --------------------------------

nuclei_qc <- function(sample_name, 
                                slide_file, 
                                nuclei_file, 
                                pdf_file,
                                force_filter = T) {
  set.seed(17)
  
  print("Reading data")
  
  print(sample_name)
  
  # We load all spots even the ones that are not over tissue
  sample_seurat <- Seurat::Load10X_Spatial(data.dir = slide_file,
                                           filename = "raw_feature_bc_matrix.h5",
                                           filter.matrix = F)
  
  nuclei_data <- read.table(file = nuclei_file,header = T, sep = ",") %>%
    rename(nuclei_count = count) %>%
    dplyr::select(barcode, nuclei_count) %>%
    column_to_rownames("barcode")
  
  sample_seurat$nuclei_count <-  nuclei_data[colnames(sample_seurat), ] 
  
  coldata <- GetTissueCoordinates(sample_seurat,
                                  cols = c("row", "col", "tissue"),
                                  scale = NULL)
  
  sample_seurat$tissue <- coldata[rownames(sample_seurat@meta.data), "tissue"]
  
  sample_seurat$tissue <- ifelse(sample_seurat$tissue == 1, 
                                 "on_tissue", 
                                 "not_on_tissue")
  
  # We do a comparison of profiles between spots in tissue and not in tissue
  
  rm(coldata)
  rm(nuclei_data)
  
  qc_a <- sample_seurat@meta.data %>%
    select(-orig.ident) %>%
    ggplot(aes(x = tissue, fill = tissue, y = nuclei_count)) +
    geom_boxplot() +
    theme_classic() +
    theme(legend.position = "none")
  
  qc_b <- sample_seurat@meta.data %>%
    ggplot(aes(y = nCount_Spatial, x = as.factor(nuclei_count))) +
    geom_boxplot() +
    theme_classic() +
    xlab("number of nuclei")
  
  qc_c <- sample_seurat@meta.data %>%
    ggplot(aes(y = nFeature_Spatial, x = as.factor(nuclei_count))) +
    geom_boxplot() +
    theme_classic() +
    xlab("number of nuclei")
  
  # Filter useful spots 
  
  sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")
  
  # Continue with analysis
  
  sample_seurat[["orig.ident"]] <- sample_name
  
  # Get mitochondrial genes percentage ------------------------------------------------
  sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, 
                                                        pattern = "^MT-")

  qc_d <- sample_seurat@meta.data %>%
    ggplot(aes(y = percent.mt, x = as.factor(nuclei_count))) +
    geom_boxplot() +
    theme_classic() +
    xlab("number of nuclei")
  
  
  # Process the data --------------------------------------------------------------------
  
  # Exclude MT and ribosomal genes ------------------------------------------------------
  
  if(force_filter) {
    
    #Ribosomal and mitochondrial genes are taken out
    mt_genes <- row.names(sample_seurat)[grepl("^MT-", row.names(sample_seurat))]
    rps_genes <- row.names(sample_seurat)[grepl("^RPS", row.names(sample_seurat))]
    mrp_genes <- row.names(sample_seurat)[grepl("^MRP", row.names(sample_seurat))]
    rpl_genes <- row.names(sample_seurat)[grepl("^RPL", row.names(sample_seurat))]
    rb_genes <- c(rps_genes, mrp_genes, rpl_genes)
    
    sample_seurat <- sample_seurat[!rownames(sample_seurat) %in% c(rb_genes, mt_genes), ]
    
    #Genes expressed in less that 10 spots are filtered
    sample_seurat <- sample_seurat[rowSums(GetAssayData(sample_seurat, assay = "Spatial") > 0) > 10, ]
    
    # Then I recalculate the number of genes and reads per spot
    sample_seurat$nFeature_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial") > 0)
    sample_seurat$nCount_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial"))
    
    filtered_spots <- WhichCells(sample_seurat,
                                 expression = nFeature_Spatial_filt > 300 & nCount_Spatial_filt > 500)
    
    sample_seurat$used_spot <- ifelse(colnames(sample_seurat) %in% filtered_spots,
                                      TRUE, FALSE)
    
    qc_e <- sample_seurat@meta.data %>%
      ggplot(aes(x = used_spot, y = nuclei_count)) +
      geom_boxplot() +
      theme_classic() +
      xlab("analyzed spots")
    
    sample_seurat <- subset(sample_seurat, 
                            subset = nFeature_Spatial_filt > 300 & 
                              nCount_Spatial_filt > 500)
    
    qc_f <- sample_seurat@meta.data %>%
      ggplot(aes(as.factor(nuclei_count))) + 
      geom_histogram(stat = "count") +
      theme_classic()
    
    qc_panel_filt <- plot_grid(qc_a, qc_b, qc_c,
                               qc_d, qc_e, qc_f, nrow = 2)
  }
  
  pdf(file = pdf_file, height = 9, width = 12.5)
  
  plot(qc_panel_filt)
  
  dev.off()
  
}

# Main -----------------------------

pwalk(param_df, nuclei_qc)





