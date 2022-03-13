# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

library(Seurat)
library(tidyverse)

sel_niches <- c("niche_0", 
                "niche_1", 
                "niche_3")

pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

niche_props <- read_csv("./results/niche_mapping/Spatial_snn_res.0.2/niche_props.csv") %>%
  left_join(pat_anns) %>%
  dplyr::filter(mol_niche %in% sel_niches,
                patient_group == "group_1") %>%
  arrange(mol_niche, -niche_prop)
  
# Plot examples

slide_name <- c("Visium_6_CK284", "AKK003_157777", "AKK002_157781", "Visium_8_CK286")
features = c("MYBPC3","ANKRD2")

for(s in slide_name) {
  
  print(s) 
  
  visium_slide <- readRDS(paste0("./processed_visium/objects/", s, ".rds"))
  
  DefaultAssay(visium_slide) <- "SCT"
  
  for(f in features) {
    
    pdf(paste0("./results/niche_mapping/Spatial_snn_res.0.2/myogenic_gene_plots/", s,"_",f, ".pdf"), height = 4.5, width = 4)
    
    plot(SpatialFeaturePlot(visium_slide, features = f))
    
    dev.off()
  }
  
  
}
