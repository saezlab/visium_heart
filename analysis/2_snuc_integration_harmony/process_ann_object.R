# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' After filtering and annotation
#' This script re-integrates data and recalculates
#' the atlas information
  
library(optparse)
library(tidyverse)
library(Seurat)
library(harmony)
library(cluster)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds object"),
  make_option(c("--out_fig_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the pdf object")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Read object to subset -----------------------------------------------------
integrated_data <- readRDS(data_path)

# Here we need to subset
integrated_data <- SplitObject(integrated_data, split.by = "orig.ident")

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
# since cell-types and niches should be greater than the number of batches

hvg_list <- map(integrated_data, function(x) {
  
  DefaultAssay(x) <- "RNA"
  
  x <- FindVariableFeatures(x, selection.method = "vst", 
                       nfeatures = 3000)
  
  x@assays[["RNA"]]@var.features
  
}) %>% unlist()

hvg_list <- table(hvg_list) %>%
  sort(decreasing = TRUE)

gene_selection <- hvg_list[1:3000] %>% names()

# We collapse again

integrated_data <- reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

# Quickly get characteristic profile of the object
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = gene_selection, 
         npcs = 30, 
         verbose = FALSE) 

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              c("orig.ident", "batch", "patient"), 
                              plot_convergence = TRUE,
                              assay.use = "RNA",
                              max.iter.harmony = 20)

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,
          reduction.name = "umap_harmony")

# Save object
saveRDS(integrated_data, file = out_path)

# Print UMAPs
umap_corrected_sample <- DimPlot(object = integrated_data, 
                                 reduction = "umap_harmony", 
                                 pt.size = .1, 
                                 group.by = "orig.ident")

umap_corrected_clustering <- DimPlot(object = integrated_data, 
                                     reduction = "umap_harmony", 
                                     pt.size = .1, 
                                     group.by = "cell_type")


pdf(file = out_fig_file, height = 10, width = 12)

print(umap_corrected_sample)
print(umap_corrected_clustering)

dev.off()

# Save UMAP info
reductions_list <-  list(meta_data = integrated_data@meta.data,
                         reduction = integrated_data@reductions[["umap_harmony"]])

saveRDS(reductions_list,
        file = gsub("[.]rds", "_umap.rds", out_path))
