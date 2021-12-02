# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we test the effects of regressing out the effects of background

library(Seurat)
library(tidyverse)
library(harmony)
library(cluster)
library(optparse)


option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--obj_out"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file")
  )

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Generate the outfiles 
pdf_out <- gsub("[.]rds", ".pdf", obj_out)

# Read the original object

integrated_data <- readRDS(data_path) 

# Get background genes from clusters in original unnannotated map

background_genes <- readRDS("./processed_snrnaseq/integration/integrated_rnasamples_mrkrs.rds")[["RNA"]] %>%
  dplyr::filter(avg_log2FC > 1, 
                p_val_adj < .000001) %>%
  arrange(cluster, -avg_log2FC) %>%
  group_by(cluster) %>%
  slice(1:30) %>%
  pull(gene) %>% unique() %>%
  sort()

# Filter gene selection to not include background genes?
# This gene selection allows not look at batch effects

gene_selection <- rownames(integrated_data@reductions$pca@feature.loadings)

gene_selection <- gene_selection[!gene_selection %in% background_genes]

background_genes <- list("background" = background_genes)

# Add module score
integrated_data <- AddModuleScore(integrated_data,
                                  assay = "RNA",
                                  features = background_genes,
                                  name = "background")

# Check, before and after integration
bckgrnd_original_umap <- FeaturePlot(integrated_data, features = "background1")
bckgrnd_original_vlns <- VlnPlot(integrated_data, features = "background1", group.by = "opt_state")
original_umap <- DimPlot(integrated_data, group.by = "opt_state")

# Copy this 
integrated_data$old_opt_state <- integrated_data$opt_state

# Re-do integration and clustering
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE, 
            vars.to.regress = "background1") %>% 
  RunPCA(features = gene_selection, 
         npcs = 30, 
         verbose = FALSE)

# Original PCs
original_pca_plt <- DimPlot(object = integrated_data, 
                            reduction = "pca",
                            pt.size = .1,
                            group.by = "orig.ident")

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              c("orig.ident", "patient", "batch"), 
                              plot_convergence = TRUE,
                              max.iter.harmony = 20)

# Corrected dimensions -----------------------
corrected_pca_plt <- DimPlot(object = integrated_data, 
                             reduction = "harmony", 
                             pt.size = .1, 
                             group.by = "orig.ident")

# Create the UMAP with new reduction -----------
integrated_data <- integrated_data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)


# Clustering and optimization -------------------------
print("Optimizing clustering")

integrated_data <- FindNeighbors(integrated_data, 
                                 reduction = "harmony", 
                                 dims = 1:30)

# Increase resolution for the lols ---------------------
seq_res <- seq(0.5, 1, 0.2)

# Delete previous clustering
integrated_data@meta.data <- integrated_data@meta.data[, !grepl("RNA_snn_res",
                                                                colnames(integrated_data@meta.data))]

# Create new clustering
integrated_data <- FindClusters(integrated_data,
                                resolution = seq_res,
                                verbose = F)


# Optimize clustering ------------------------------------------------------
cell_dists <- dist(integrated_data@reductions$harmony@cell.embeddings,
                   method = "euclidean")

cluster_info <- integrated_data@meta.data[,grepl("RNA_snn_res",
                                                 colnames(integrated_data@meta.data))] %>%
  dplyr::mutate_all(as.character) %>%
  dplyr::mutate_all(as.numeric)

silhouette_res <- apply(cluster_info, 2, function(x){
  si <- silhouette(x, cell_dists)
  mean(si[, 'sil_width'])
})

integrated_data[["opt_state"]] <- integrated_data[[names(which.max(silhouette_res))]]

Idents(integrated_data) = "opt_state"

# Final plots

pdf(pdf_out)

# Plot OLD UMAPs and clustering
plot(original_umap)
plot(bckgrnd_original_umap)
plot(bckgrnd_original_vlns)

# Plot: If this looks the same then this method wasn't useful
umap_corrected_clustering <- DimPlot(object = integrated_data, 
                                     reduction = "umap", 
                                     pt.size = .1, 
                                     group.by = "opt_state")

plot(umap_corrected_clustering)

umap_corrected_oldclustering <- DimPlot(object = integrated_data, 
                                        reduction = "umap", 
                                        pt.size = .1, 
                                        group.by = "old_opt_state")

plot(umap_corrected_oldclustering)

plot(FeaturePlot(integrated_data, features = "background1"))
plot(VlnPlot(integrated_data, features = "background1", group.by = "opt_state"))

dev.off()

# Save object ------------------------------------------------------

saveRDS(integrated_data, file = obj_out)