library(tidyverse)
library(Seurat)


ctrls <- read_csv("./markers/visium_patient_anns_revisions.csv") %>%
  dplyr::filter(major_labl == "CTRL") %>%
  pull(sample_id)

walk(ctrls, function(slide_id) {
  
  slide_file <- paste0("./processed_visium/objects/", slide_id, ".rds")
  pdf_file <- paste0("./results/niche_mapping/Spatial_snn_res.0.2/myogenic_gene_plots/", slide_id)
  
  visium_slide <- readRDS(slide_file)
  
  DefaultAssay(visium_slide) <- "SCT"
  
  walk(c("MYBPC3", "ANKRD2", "RYR2"), function(g){
    
    plt <- SpatialFeaturePlot(visium_slide, features = g,max.cutoff = "q99", min.cutoff = "q1")
    
    pdf(paste0(pdf_file,"_CTRL", "_", g, ".pdf"), height = 4, width = 4)
    
    plot(plt)
    
    dev.off()
    
  })
  
  
})

bzs <- read_csv("./markers/visium_patient_anns_revisions.csv") %>%
  dplyr::filter(major_labl == "BZ") %>%
  pull(sample_id)

walk(bzs, function(slide_id) {
  
  slide_file <- paste0("./processed_visium/objects/", slide_id, ".rds")
  pdf_file <- paste0("./results/niche_mapping/Spatial_snn_res.0.2/myogenic_gene_plots/", slide_id)
  
  visium_slide <- readRDS(slide_file)
  
  DefaultAssay(visium_slide) <- "SCT"
  
  walk(c("MYBPC3", "ANKRD2", "RYR2"), function(g){
    
    plt <- SpatialFeaturePlot(visium_slide, features = g,max.cutoff = "q99", min.cutoff = "q1")
    
    pdf(paste0(pdf_file,"_bz", "_", g, ".pdf"), height = 4, width = 4)
    
    plot(plt)
    
    dev.off()
    
  })
  
  
})

visium_slide <- readRDS("./processed_visium/objects/AKK002_157781.rds")  %>%
  positive_states(., assay = state_origin) %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.1)



SpatialFeaturePlot(visium_slide,c("CCDC80", "TMSB10"))



