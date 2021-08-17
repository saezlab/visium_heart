# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot the niche characteristics

library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(corrplot)

degs <- readRDS("./processed_visium/integration/integrated_slides_dea.rds")

visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate(niche_ct_corr = map(visium_file, function(f_path) {
    
    visium_slide <- readRDS(f_path)
    Idents(visium_slide) <- "opt_clust_integrated"

    niche_cors <- map(set_names(levels(Idents(visium_slide))),
                       function(niche) {
                         
                         niche_cells <- WhichCells(visium_slide, idents = c(niche))
                         
                         prop_mat <- GetAssayData(visium_slide, 
                                                  assay = "c2l_major_props") %>%
                           as.matrix()
                         
                         prop_mat <- prop_mat[, niche_cells]
                         
                         cor_mat <- t(prop_mat) %>%
                           cor(method = "spearman")
                       })
    
  })) %>%
  unnest_longer("niche_ct_corr")


visium_df <- visium_df  %>%
  mutate(niche_ct_corr = map(niche_ct_corr, function(x) {
    
    x[upper.tri(x, diag = T)] <- NA 
    
    x %>%
      as.data.frame() %>%
      rownames_to_column("cellA") %>%
      pivot_longer(-cellA, 
                   values_to = "spearman_cor",
                   names_to = "cellB") %>%
      na.omit
    
    
  }))


allcors_long <- visium_df %>%
  unnest(niche_ct_corr) %>%
  ungroup() %>%
  arrange(niche_ct_corr_id, cellA, cellB) %>%
  mutate(cor_id = paste0(cellA, "_", cellB),
         niche_id = paste0(niche_ct_corr_id, "_", sample)) %>%
  dplyr::select(spearman_cor, cor_id, niche_id) %>%
  ungroup()

allcors_mat <- allcors_long %>%
  pivot_wider(values_from = spearman_cor, names_from = cor_id) %>%
  column_to_rownames("niche_id") %>%
  as.matrix()

order_rows <- hclust(dist(allcors_mat))$order
order_rows_labs <- rownames(allcors_mat)[order_rows]
order_cols <- hclust(dist(t(allcors_mat)))$order
order_cols_labs <- colnames(allcors_mat)[order_cols]

ggplot(allcors_long, 
       aes(y = factor(cor_id, levels = order_cols_labs),
           x = factor(niche_id, levels = order_rows_labs),
           fill = spearman_cor)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))

Heatmap(scale(t(allcors_mat)))

###############

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate(niche_ct_meanprops = map(visium_file, function(f_path) {
    
    visium_slide <- readRDS(f_path)
    Idents(visium_slide) <- "opt_clust_integrated"
    
    niche_cors <- map(set_names(levels(Idents(visium_slide))),
                      function(niche) {
                        print(f_path)
                        print(niche)
                        
                        niche_cells <- WhichCells(visium_slide, idents = c(niche))
                        
                        prop_mat <- GetAssayData(visium_slide, 
                                                 assay = "c2l_major_props") %>%
                          as.matrix()
                        
                        print(length(niche_cells))
                        
                        prop_mat <- prop_mat[, niche_cells, drop = F] %>%
                          rowMeans()
                      })
    
  })) %>%
  unnest_longer("niche_ct_meanprops")

visium_df <- visium_df %>%
  dplyr::select(sample, niche_ct_meanprops, niche_ct_meanprops_id) %>%
  mutate(niche_ct_meanprops = map(niche_ct_meanprops, function(x) {
    
    df_props <- x %>%
      as.data.frame()
    
    colnames(df_props) = "mean_prop"
    
    df_props %>% rownames_to_column("cell_type")
    
  })) %>% unnest() %>%
  arrange(niche_ct_meanprops_id, cell_type, sample) %>%
  mutate(niche_id = paste0(sample, "_", niche_ct_meanprops_id))

props_mat <- visium_df %>%
  dplyr::select(cell_type, niche_id, mean_prop) %>%
  pivot_wider(values_from = mean_prop,
              names_from = cell_type) %>%
  column_to_rownames("niche_id") %>%
  as.matrix()

ggplot(visium_df, aes(x = niche_id, 
                      y = cell_type, 
                      fill = mean_prop)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

pdf("./processed_visium/integration/niche_compositions/mean_props.pdf",
    width = 8, height = 16)

draw(Heatmap(scale(props_mat), cluster_rows = F,
             row_split = factor(rownames(props_mat) %>% strsplit("[1-9]_") %>% map_chr(., ~.x[2]))))

dev.off()
