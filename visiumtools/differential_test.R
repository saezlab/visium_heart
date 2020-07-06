# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#'Performs wilcoxon test to identify differential features
#'of identities (by default from clustering) of selected slides

find_allfeat = function(visium_slide, assays_collection = c("SCT","dorothea","progeny",
                                                            "ctscores")){
  require(Seurat)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(dplyr)
  
  possible_assays = set_names(assays_collection)
  
  diff_features = lapply(possible_assays, function(x){
    
    DefaultAssay(visium_slide) = x
    
    vis_wilcox_markers = FindAllMarkers(visium_slide, 
                                        logfc.threshold = 0.05,
                                        test.use = "wilcox")
  })
  
  return(diff_features)
  
} 



