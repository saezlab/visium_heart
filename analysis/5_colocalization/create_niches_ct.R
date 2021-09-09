# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We integrate the deconvolution results of all slides and create a UMAP
#' representation of this,
#' 
#' We compare ILR transformations of the data and pure proportions in the definition of the UMAP

library(scater)
library(tidyverse)
library(uwot)
library(cowplot)
library(factoextra)
library(Seurat)
library(ggrastr)
library(compositions)

c2l_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"
assay_name <- "c2l"
q05 <- TRUE
make_proportions <- TRUE

# Get patient meta-data
sample_dict <- read.table("./markers/visium_annotations_ext.txt", 
                          sep = "\t", header = T)

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

# Plot UMAP
comp_umap <- umap(integrated_compositions, 
                 n_neighbors = 30, n_epochs = 1000, 
                 metric = "cosine") %>%
  as.data.frame() %>%
  mutate(row_id = rownames(integrated_compositions))

comp_umap %>%
  left_join(atlas_meta, by = c("row_id")) %>%
  ggplot(aes(x = V1, y = V2, 
             color = opt_clust_integrated)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2")

# Plot UMAP with cell_scores

cts <- set_names(colnames(integrated_compositions))

pdf("./results/niche_mapping/ct_niches/ct_proportions_umap.pdf", height = 6, width = 7)

walk(cts, function(ct){
  
  plot_df <- comp_umap %>%
    mutate(ct_prop = integrated_compositions[ , ct])
  
  plt <- plot_df %>%
    ggplot(aes(x = V1, y = V2, 
               color = ct_prop)) +
    ggrastr::geom_point_rast(size = 0.07) +
    theme_classic() +
    ggtitle(ct) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_colour_gradient(low = "#ffd89b", high = "#19547b", limits = c(0,1))
  
  plot(plt)
  
})

dev.off()

# Create UMAP
comp_umap <- umap(integrated_compositions, 
                  n_neighbors = 30, n_epochs = 1000, 
                  metric = "cosine") %>%
  as.data.frame() %>%
  mutate(row_id = rownames(integrated_compositions))

comp_umap %>%
  left_join(atlas_meta, by = c("row_id")) %>%
  ggplot(aes(x = V1, y = V2, 
             color = opt_clust_integrated)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2")

# Plot UMAP with cell_scores

cts <- set_names(colnames(integrated_compositions))

pdf("./results/niche_mapping/ct_niches/ct_proportions_umap.pdf", height = 6, width = 7)

walk(cts, function(ct){
  
  plot_df <- comp_umap %>%
    mutate(ct_prop = integrated_compositions[ , ct])
  
  plt <- plot_df %>%
    ggplot(aes(x = V1, y = V2, 
               color = ct_prop)) +
    ggrastr::geom_point_rast(size = 0.07) +
    theme_classic() +
    ggtitle(ct) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_colour_gradient(low = "#ffd89b", high = "#19547b", limits = c(0,1))
  
  plot(plt)
  
})

dev.off()

# Create PCA

comp_pcs_res <- prcomp(x = integrated_compositions)
comp_pcs <- comp_pcs_res$x %>%
  as.data.frame() %>%
  rownames_to_column("row_id")

comp_pcs %>%
  left_join(atlas_meta, by = "row_id") %>%
  ggplot(aes(x = PC1, y = PC2, 
             color = opt_clust_integrated)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("PC_1") +
  ylab("PC_2")

# Compositions shouldn't be treated like this
# Let's use isometric log ratios

baseILR <- ilrBase(x = integrated_compositions,
                   method = "basic")

cell_ilr <- as.matrix(ilr(integrated_compositions, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

comp_pcs_res <- prcomp(x = cell_ilr)
comp_pcs <- comp_pcs_res$x %>%
  as.data.frame() %>%
  rownames_to_column("row_id")

comp_pcs %>%
  left_join(atlas_meta, by = "row_id") %>%
  ggplot(aes(x = PC1, y = PC2, 
             color = opt_clust_integrated)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("PC_1") +
  ylab("PC_2")


# What about UMAP on this?
comp_umap <- umap(cell_ilr, 
                  n_neighbors = 30, n_epochs = 1000) %>%
  as.data.frame() %>%
  mutate(row_id = rownames(cell_ilr))

comp_umap %>%
  left_join(atlas_meta, by = c("row_id")) %>%
  ggplot(aes(x = V1, y = V2, 
             color = opt_clust_integrated)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2")

pdf("./results/niche_mapping/ct_niches/ct_ILR_umap.pdf", height = 6, width = 7)

walk(cts, function(ct){
  
  plot_df <- comp_umap %>%
    mutate(ct_prop = integrated_compositions[ , ct])
  
  plt <- plot_df %>%
    ggplot(aes(x = V1, y = V2, 
               color = ct_prop)) +
    ggrastr::geom_point_rast(size = 0.07) +
    theme_classic() +
    ggtitle(ct) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_colour_gradient(low = "#ffd89b", high = "#19547b", limits = c(0,1))
  
  plot(plt)
  
})

dev.off()

# What about clustering?

saveRDS(list("ilr_counts" = cell_ilr,
             "meta_data" = atlas_meta), 
        file = "./processed_visium/integration/integrated_slides_ct_ILR.rds")



