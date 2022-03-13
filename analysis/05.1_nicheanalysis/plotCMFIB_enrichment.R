# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Plot only fibroblasts and CM states vs niche
library(tidyverse)
library(ComplexHeatmap)

make_hmaps <- function(table_file, res_dir, alias) {
  
  stats_table <- read_csv(table_file) 
  
  # Fib
  
  enrich_lab <- stats_table %>%
    dplyr::filter(grepl("Fib", name)) %>% 
    dplyr::select(name, mol_niche, sign) %>%
    dplyr::mutate(sign = ifelse(is.na(sign), "", sign)) %>%
    pivot_wider(names_from = mol_niche, values_from = sign) %>%
    column_to_rownames("name") %>%
    as.matrix() 
  
  enrich_val <- stats_table %>%
    dplyr::filter(grepl("Fib", name)) %>%
    group_by(name) %>%
    #mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop)) %>%
    dplyr::select(name, mol_niche, scaled_median_prop) %>%
    pivot_wider(names_from = mol_niche, values_from = scaled_median_prop) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  fib_plt <- Heatmap(enrich_val, name = "scaled comp", rect_gp = gpar(col = "black", lwd = 1),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf(enrich_lab[i, j]), x, y, gp = gpar(fontsize = 10))
                     })
  
  pdf(paste0(res_dir, "/", "Fib", alias,".pdf"), height = 2.5, width = 6)
  
  ComplexHeatmap::draw(fib_plt)
  
  dev.off()
  
  # CM
  
  enrich_lab <- stats_table %>%
    dplyr::filter(grepl("CM", name)) %>% 
    dplyr::select(name, mol_niche, sign) %>%
    dplyr::mutate(sign = ifelse(is.na(sign), "", sign)) %>%
    pivot_wider(names_from = mol_niche, values_from = sign) %>%
    column_to_rownames("name") %>%
    as.matrix() 
  
  enrich_val <- stats_table %>%
    dplyr::filter(grepl("CM", name)) %>%
    group_by(name) %>%
    mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop)) %>%
    dplyr::select(name, mol_niche, scaled_median_prop) %>%
    pivot_wider(names_from = mol_niche, values_from = scaled_median_prop) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  cm_plt <- Heatmap(enrich_val, name = "scaled comp", rect_gp = gpar(col = "black", lwd = 1),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf(enrich_lab[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
  
  
  
  pdf(paste0(res_dir, "/", "CM", alias,".pdf"), height = 2.5, width = 6)
  
  ComplexHeatmap::draw(cm_plt)
  
  dev.off()
  
  #Myeloid
  
  enrich_lab <- stats_table %>%
    dplyr::filter(grepl("Myeloid", name)) %>% 
    dplyr::select(name, mol_niche, sign) %>%
    dplyr::mutate(sign = ifelse(is.na(sign), "", sign)) %>%
    pivot_wider(names_from = mol_niche, values_from = sign) %>%
    column_to_rownames("name") %>%
    as.matrix() 
  
  enrich_val <- stats_table %>%
    dplyr::filter(grepl("Myeloid", name)) %>%
    group_by(name) %>%
    mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop)) %>%
    dplyr::select(name, mol_niche, scaled_median_prop) %>%
    pivot_wider(names_from = mol_niche, values_from = scaled_median_prop) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  cm_plt <- Heatmap(enrich_val, name = "scaled comp", rect_gp = gpar(col = "black", lwd = 1),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf(enrich_lab[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
  
  
  
  pdf(paste0(res_dir, "/", "Myeloid", alias,".pdf"), height = 2.5, width = 6)
  
  ComplexHeatmap::draw(cm_plt)
  
  dev.off()
}

# Cell-type

int_res <- c("composition_niche",
             "Spatial_snn_res.0.2")


make_hmaps(table_file = "./results/niche_mapping/composition_niche/cell_state_enrich_wilcox.csv",
           res_dir = "./results/niche_mapping/composition_niche/",
           alias = "composition_niche")

make_hmaps(table_file = "./results/niche_mapping/Spatial_snn_res.0.2/cell_state_enrich_wilcox.csv",
           res_dir = "./results/niche_mapping/Spatial_snn_res.0.2/",
           alias = "Spatial_snn_res.0.2")



