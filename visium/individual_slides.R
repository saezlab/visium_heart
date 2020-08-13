# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing visium slides
#' Adds PROGENy scores, TF activities,
#' Independent cell type markers

source("./visium_exploratory/slide_processing.R")

#Markers from single cell - FROM VISIUM dataset
scell_marker_files = list.files("./results/single_sample/all_markers")
scell_marker_files = unlist(lapply(strsplit(scell_marker_files,"\\."), 
                                   function(x) x[1]))

sample_dictionary = readRDS("./sample_dictionary.rds")

slides_ids = c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785")

for(slide in slides_ids){
  print(slide)
  
  out_path = sprintf("mkdir ./results/single_slide/%s",slide)
  system(out_path)
  
  #Processing visium slide
  slide_seurat = process_visium(dir_path = sprintf("./visium_data/%s", slide))
  
  
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
  
  #ECM scores
  DefaultAssay(slide_seurat) = "ECM"
  plot(SpatialFeaturePlot(object = slide_seurat,
                          features = rownames(slide_seurat)))
  
  
  dev.off()
  
  saveRDS(slide_seurat, file = sprintf("./results/single_slide/%s/%s.rds",
                                       slide,slide))
  
}

## Add matched cell-type scores as an extra step
for(slide in slides_ids){
  print(slide)
  # File definition
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  #Processing visium slide
  slide_seurat = readRDS(slide_file)
  
  #Cell type scores from specific scRNA if data available
  scell_id = dplyr::filter(sample_dictionary, visium_ids == slide) %>%
    dplyr::select(scell_ids) %>% pull()
  
  if(scell_id %in% scell_marker_files){
    
    # Here create score matrix
    scell_markers = read.table(file = sprintf("./results/single_sample/all_markers/%s.AllMarkers.txt",
                                              scell_id),
                               header = T, sep = "\t",
                               stringsAsFactors = F) %>%
      dplyr::mutate(cluster = gsub(" ","_",cluster)) %>%
      dplyr::mutate(cluster = gsub("\\+","pos",cluster)) %>%
      group_by(cluster) %>% arrange(p_val_adj) %>%
      slice(1:100) %>% dplyr::select(avg_logFC,cluster,gene) %>%
      pivot_wider(names_from = cluster,values_from = avg_logFC)
    
    markers_mat = as.matrix(scell_markers[,-1])
    rownames(markers_mat) = scell_markers$gene
    markers_mat[is.na(markers_mat)] = 0
    
    exprdata = as.matrix(slide_seurat[["SCT"]]@data)
    genes_in_markers = rownames(exprdata) %in% rownames(markers_mat)
    exprdata = exprdata[genes_in_markers,]
    markers_mat = markers_mat[rownames(exprdata),]
    
    ct_scores = scale(t(t(markers_mat) %*% exprdata))
    
    slide_seurat[['ctscores_match']] = CreateAssayObject(counts = t(ct_scores))
    
    DefaultAssay(slide_seurat) = "ctscores_match"
    
    saveRDS(slide_seurat, file = sprintf("./results/single_slide/%s/%s.rds",
                                         slide,slide))
    
  }
  
}
