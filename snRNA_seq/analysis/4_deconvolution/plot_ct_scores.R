library(Seurat)
library(tidyverse)
library(cowplot)


visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

assays <- set_names(c("c2l_major", "c2l_states", 
            "c2l_major_props", "c2l_states_props"))

SpatialPal = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))


visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate(assay_plots = map(visium_file, function(f_path) {
    
    visium_slide <- readRDS(f_path)
    
    assay_plts <- map(assays, function(assay_name) {
      
      features <- rownames(GetAssayData(visium_slide,
                                        slot = "data",
                                        assay = assay_name))
      
      DefaultAssay(visium_slide) <- assay_name
      
      feat_plots <- map(features, function(spec_f) {
        
        f_p <- SpatialFeaturePlot(visium_slide, 
                           features = spec_f)
        
        if(grepl(pattern = "props", assay_name)) {
          f_p <- f_p + 
            scale_fill_gradientn(colours = SpatialPal(length(seq(0,1,.1))),
                                 limits = c(0,1),
                                 breaks = c(0, 0.25, 0.5, 0.75, 1))
        }
        
        return(f_p)
        
      })
      
      niche_plots <- plot_grid(plotlist= feat_plots, ncol = 4)
      
    })
    
  }))


visium_df <- visium_df %>%
  unnest_longer(assay_plots) %>%
  group_by(assay_plots_id) %>%
  nest() %>%
  dplyr::mutate(out_pdf = paste0("./processed_visium/ct_plots/", assay_plots_id, ".pdf"))

walk2(visium_df$out_pdf, visium_df$data, function(x, y) {
  
  pdf(file = x, height = 15, width = 15)
  
  walk(y, print)
  
  dev.off()
})


