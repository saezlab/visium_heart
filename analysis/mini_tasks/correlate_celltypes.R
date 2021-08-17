library(Seurat)
library(tidyverse)
library(cowplot)
library(ComplexHeatmap)


visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

assays <- set_names(c("c2l_major", "c2l_states", 
                      "c2l_major_props", "c2l_states_props"))

SpatialPal = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))


visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate(assay_cormats = map(visium_file, function(f_path) {
    
    visium_slide <- readRDS(f_path)
    
    assay_cors<- map(assays, function(assay_name) {
      
      cor_mat <- as.matrix(GetAssayData(visium_slide,
                                        slot = "data",
                                        assay = assay_name)) %>%
        t() %>%
        cor(method = "spearman")
      
      
    })
    
  }))

visium_df <- visium_df %>% 
  unnest_longer("assay_cormats") %>%
  group_by(assay_cormats_id) %>%
  nest() %>%
  dplyr::mutate(out_pdf = paste0("./processed_visium/ct_plots/cormat", 
                                 assay_cormats_id, 
                                 ".pdf"))

walk2(visium_df$out_pdf, visium_df$data, function(x, y) {
  
  pdf(file = x, height = 15, width = 15)
  
  walk2(y$assay_cormats, y$sample, function(assay_cormats, sample){
    
    hmap_plt <- ComplexHeatmap::Heatmap(assay_cormats,
                                        name = sample,
                                        show_row_dend = F)
    draw(hmap_plt)
    
  })
  
  dev.off()
})





