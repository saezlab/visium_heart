# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing visium slides
#' 
#' 

source("./visium_exploratory/slide_processing.R")

slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

for(slide in slides_ids){
  print(slide)
  
  out_path = sprintf("mkdir ./results/single_slide/%s",slide)
  system(out_path)
  
  #Processing visium slide
  slide_seurat = process_visium(dir_path = sprintf("./visium_data/%s", slide))
  saveRDS(slide_seurat, file = sprintf("./results/single_slide/%s/%s.rds",
                                       slide,slide))
  
  #Creating output pdf
  pdf(file = sprintf("./results/single_slide/%s/feature_plot.pdf",
                     slide), onefile=TRUE,width = 11, height = 10)
  
  #Counts
  SpatialFeaturePlot(slide_seurat, features = "nCount_Spatial") + 
    theme(legend.position = "right")
  
  
  #Cell type scores
  DefaultAssay(slide_seurat) = "ctscores"
  plot(SpatialFeaturePlot(object = slide_seurat,
                     features = rownames(slide_seurat)))
  
  #PROGENy scores
  DefaultAssay(slide_seurat) = "progeny"
  plot(SpatialFeaturePlot(object = slide_seurat,
                          features = rownames(slide_seurat)))
  
  dev.off()
  
  
}


