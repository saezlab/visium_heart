# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate cell2location results from all slides
#' 
library(tidyverse)
library(ComplexHeatmap)

c2l_folder <- "./results/deconvolution_models/location_models/density_tables_rds/"

# Get cell2location files --------------------------------
c2l_files <- list.files(c2l_folder, full.names = F)

c2l_samples <- map_chr(strsplit(c2l_files,".rds"), 
                       ~ .x[1])

c2l_df <- tibble(c2l_file = paste0(c2l_folder, 
                                   c2l_files),
                 sample = c2l_samples) 

# Generates list of matrix of c2l results
list_matrices <- map2(c2l_df$c2l_file, c2l_df$sample, function(f, s) {
  mat <- readRDS(f)
  rownames(mat) <- paste0(s, "..", rownames(mat))
  prop_mat <- base::apply(mat, 1, function(x) {
    
    x/sum(x)
    
  })
  
  return(t(prop_mat))
})

integrated_compositions <- reduce(list_matrices, rbind) 

cor_mat <- cor(integrated_compositions)
cor_mat_order <- hclust(as.dist(1-cor_mat))
cor_mat_order <- cor_mat_order$labels[cor_mat_order$order]

cor_mat <- cor_mat[cor_mat_order,cor_mat_order]

cor_mat[lower.tri(cor_mat,diag = T)] <- NA

cor_plt_dat <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column("cell_a") %>%
  pivot_longer(-cell_a, values_to = "p_corr", names_to = "cell_b") %>%
  na.omit() %>%
  dplyr::mutate(cell_a = factor(cell_a,
                                levels = cor_mat_order),
                cell_b = factor(cell_b,
                                levels = cor_mat_order))


cor_plt <- ggplot(cor_plt_dat, aes(x = cell_a, y = cell_b, fill = p_corr)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
        axis.text = element_text(size = 12)) +
  scale_fill_gradient2() +
  coord_equal()

pdf("./results/tissue_structure/colocalization/c2l_correlation.pdf", height = 4, width = 5)

plot(cor_plt)

dev.off()
