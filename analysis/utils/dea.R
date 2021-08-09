# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

library(Seurat)
library(tidyverse)

#'Performs wilcoxon test to identify differential features
#'of identities (by default from clustering) of selected assays
#' @param visium_slide: Seurat object with defined identities
#' @param assay_collection: vector with names of assays to test
#' @param logfc.threshold: LFC threshold to use
#' @return a list with the results for each assay defined
find_allfeat <- function(visium_slide, 
                        assays_collection = c("dorothea",
                                              "progeny",
                                              "RNA"),
                        logfc.threshold = 0.05,
                        only.pos = FALSE){
  
  possible_assays <- set_names(assays_collection)
  
  diff_features <- map(possible_assays, function(x) {
    
    DefaultAssay(visium_slide) = x
    
    vis_wilcox_markers = FindAllMarkers(visium_slide, 
                                        logfc.threshold = logfc.threshold,
                                        test.use = "wilcox",
                                        only.pos = only.pos)
  })
  
  return(diff_features)
} 

#'Performs a hypergeometric test to enrich a gene list into an annotation data base
#' @param geneList: a vector containing gene names
#' @param Annotation_DB: a list of gene_sets that contains a vector of genes per gene_set
#' @return a data frame with overrepresented genes

GSE_analysis = function(geneList,Annotation_DB){
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  ResultsDF = ResultsDF %>% 
    rownames_to_column("gset") %>% 
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"), 
              as.numeric) %>% 
    dplyr::arrange(corr_p_value,GenesInList)
  
  return(ResultsDF)
  
}