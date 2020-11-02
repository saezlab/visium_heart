# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing visium slides
#' Adds PROGENy scores, TF activities,
#' Independent cell type markers

source("./visium/visiumtools/reading2clustering.R")
source("./visium/visiumtools/funcomics.R")
source("./visium/visiumtools/ECM_immune_scores.R")

slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

for(slide in slides_ids){
  print(slide)
  
  out_path = sprintf("mkdir ./results/single_slide/%s",slide)
  system(out_path)
  
  #Processing visium slide - visium tools
  slide_seurat = process_visium(dir_path = sprintf("./visium_data/%s", slide),
                                var_features = "seurat", 
                                out_dir = sprintf("./results/single_slide/%s/%s_vtools_processing.pdf",slide,slide),
                                resolution = 1,verbose = FALSE)
  
  #Adding dorothea - progeny
  slide_seurat = add_funcomics(visium_slide = slide_seurat, 
                               top = 1000)
    
  #Adding ECM scores
  slide_seurat[['ECM']] = CreateAssayObject(counts = get_ECMscores(seurat_obj = slide_seurat,
                                                          NABA_SETS = NABA_SETS))
  
  #Adding 
  
  #Creating output pdf
  pdf(file = sprintf("./results/single_slide/%s/feature_plot.pdf",
                     slide), onefile=TRUE,width = 11, height = 10)
  
  #Counts
  SpatialFeaturePlot(slide_seurat, features = "nCount_Spatial") + 
    theme(legend.position = "right")
  
  
  #PROGENy scores
  DefaultAssay(slide_seurat) = "progeny"
  plot(SpatialFeaturePlot(object = slide_seurat,
                          features = rownames(slide_seurat)))
  
  #ECM scores
  DefaultAssay(slide_seurat) = "ECM"
  plot(SpatialFeaturePlot(object = slide_seurat,
                          features = rownames(slide_seurat)))
  
  dev.off()
  
  saveRDS(slide_seurat, file = sprintf("./results/single_slide/%s/%s.rds",
                                       slide,slide))
  
}
