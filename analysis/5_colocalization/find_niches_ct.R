# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we perform iterative clustering of the ILR data in a Seurat object to define niches
#' 1) Get composition matrix
#' 2) Do ILR transformation
#' 3) SNN and Louvain
#' 4) Show the cell-type specific responses per niche

library(compositions)
library(Seurat)
library(tidyverse)
library(clustree)
library(uwot)

# Here we will get the matrix of compositions using c2l results

c2l_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"
assay_name <- "c2l"

# Define resolution - You have to do this manually
niche_resolution <- "ILR_snn_res.1"

# Get patient meta-data -----------------------------------------------------------------
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

# Generate ILR transformation
baseILR <- ilrBase(x = integrated_compositions,
                   method = "basic")

cell_ilr <- as.matrix(ilr(integrated_compositions, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

# Make Seurat object
seurat_ilr <- CreateSeuratObject(t(cell_ilr %>% as.data.frame() %>% as.matrix() ), 
                                 project = "SeuratProject", 
                                 assay = "ILR",
                                 in.cells = 0, min.features = 0, names.field = 1,
                                 names.delim = "_", meta.data = atlas_meta %>% column_to_rownames("row_id"))

# Neighbor graph

seurat_ilr <- ScaleData(seurat_ilr, 
                           verbose = FALSE, 
                           features = rownames(seurat_ilr))

seurat_ilr <- RunPCA(seurat_ilr,
                     verbose = FALSE,
                     features = rownames(seurat_ilr),
                     approx=FALSE)

seurat_ilr <- FindNeighbors(seurat_ilr)

seurat_ilr@graphs$ILR_snn <- neighbor_graph$snn
seurat_ilr@graphs$ILR_nn <- neighbor_graph$nn

print("Created network")

seurat_ilr <- FindClusters(seurat_ilr, resolution = seq(0.4, 1.6, 0.2))

clustree_plt <- clustree(seurat_ilr, 
                         prefix = paste0(DefaultAssay(seurat_ilr), 
                                         "_snn_res."))

pdf("./results/niche_mapping/ct_niches/ct_ILR_clustree.pdf", height = 15, width = 15)

plot(clustree_plt)

dev.off()

saveRDS(seurat_ilr, file = "./processed_visium/integration/integrated_slides_ct_ILR_srt.rds")


# Create UMAP and plot the compositions

comp_umap <- umap(cell_ilr, 
                  n_neighbors = 30, n_epochs = 1000) %>%
  as.data.frame() %>%
  mutate(row_id = rownames(cell_ilr))

comp_umap <- comp_umap %>%
  left_join(seurat_ilr@meta.data %>%
              rownames_to_column("row_id"), by ="row_id")



pdf("./results/niche_mapping/ct_niches/ct_ILR_umap.pdf", height = 6, width = 7)

plt <- comp_umap %>%
  ggplot(aes(x = V1, y = V2, 
             color = ILR_snn_res.1)) +
  ggrastr::geom_point_rast(size = 0.1) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2")

plot(plt)

cts <- set_names(colnames(integrated_compositions))

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

cluster_cols <- grepl("snn", colnames(comp_umap))
cluster_cols <- colnames(comp_umap)[cluster_cols]
cluster_info <- comp_umap[,c("row_id", "orig.ident", cluster_cols)]

# What are the cells that define the niches?

integrated_compositions <- integrated_compositions %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id,values_to = "ct_prop", names_to = "cell_type")

integrated_compositions <- integrated_compositions %>%
  left_join(cluster_info)

# Select a niche resolution
niche_summary_pat <- integrated_compositions %>%
  dplyr::select_at(c("row_id", "cell_type", "ct_prop", "orig.ident", niche_resolution)) %>%
  rename("ct_niche" = niche_resolution) %>%
  mutate(ct_niche = paste0("niche_", ct_niche)) %>%
  group_by(orig.ident, ct_niche, cell_type) %>%
  summarize(mean_ct_prop = mean(ct_prop))
  
niche_summary <- niche_summary_pat %>%
  ungroup() %>%
  group_by(ct_niche, cell_type) %>%
  summarise(patient_mean_ct_prop = mean(mean_ct_prop))

# Data manipulation to have clustered data

niche_summary_mat <- niche_summary %>%
  pivot_wider(values_from = patient_mean_ct_prop, 
              names_from =  cell_type, values_fill = 0) %>%
  column_to_rownames("ct_niche") %>%
  as.matrix()

niche_order <- hclust(dist(niche_summary_mat))
niche_order <- niche_order$labels[niche_order$order]

ct_order <- hclust(dist(t(niche_summary_mat)))
ct_order <- ct_order$labels[ct_order$order]

# Find characteristic cell types of each niche
# We have per patient the proportion of each cell-type in each niche

run_wilcox_all <- function(prop_data) {
  
  prop_data_group <- prop_data[["ct_niche"]] %>%
    unique() %>%
    set_names()
  
  map(prop_data_group, function(g) {
    
    test_data <- prop_data %>%
      mutate(test_group = ifelse(ct_niche == g,
                                 "target", "rest"))
    
    wilcox.test(mean_ct_prop ~ test_group, 
                data = test_data,
                alternative = "two.sided") %>%
      broom::tidy()
  }) %>% enframe("ct_niche") %>%
    unnest()
  
}

wilcoxon_res <- niche_summary_pat %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = map(data, run_wilcox_all)) %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value))

wilcoxon_res <- wilcoxon_res %>%
  mutate(significant = ifelse(p_corr <= 0.15, "*", ""))

mean_ct_prop_plt <- niche_summary %>%
  left_join(wilcoxon_res, by = c("ct_niche", "cell_type")) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order),
         ct_niche = factor(ct_niche, levels = niche_order)) %>%
  ggplot(aes(x = cell_type, y = ct_niche, fill = patient_mean_ct_prop)) +
  geom_tile() +
  geom_text(aes(label = significant)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.y = element_text(size=12)) +
  scale_fill_gradient(high = "#ffd89b", low = "#19547b", limits = c(0,1)) 

# Finally describe the proportions of those niches in all the data
cluster_counts <- cluster_info %>%
  dplyr::select_at(c("row_id", "orig.ident", niche_resolution)) %>%
  rename("ct_niche" = niche_resolution) %>%
  mutate(ct_niche = paste0("niche_", ct_niche)) %>%
  group_by(ct_niche) %>%
  summarise(nspots = length(ct_niche)) %>%
  mutate(prop_spots = nspots/sum(nspots))

barplts <- cluster_counts %>%
  mutate(ct_niche = factor(ct_niche, levels = niche_order)) %>%
  ggplot(aes(y = ct_niche, x = prop_spots)) +
  geom_bar(stat = "identity") +
  theme_classic() + ylab("") +
  theme(axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(size=12))

niche_summary_plt <- cowplot::plot_grid(mean_ct_prop_plt, barplts, align = "hv", axis = "tb")

pdf("./results/niche_mapping/ct_niches/characteristic_ct_niches.pdf", height = 6, width = 6)

plot(niche_summary_plt)

dev.off()


saveRDS(cluster_info, "./results/niche_mapping/ct_niches/niche_annotation_ct.rds")
