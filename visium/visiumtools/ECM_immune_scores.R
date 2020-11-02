# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing of breast cancer visium slides
#' It normalizes the data using SCT transform,
#' calculates TF and PROGENy activities and
#' defines a matrix of expressed ligands
#' 

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(viper)
library(progeny)
library(ggplot2)
library(MISTy)
# General data

#NABA genesets
processNABA = function(filepath = "./markers/NABAgsets.xls") {
  con = file(filepath, "r")
  naba_gsets = list()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    split_line = unlist(strsplit(line,split="\t"))
    naba_gsets[[split_line[1]]] = split_line[3:length(split_line)]
  }
  close(con)
  return(naba_gsets)
}

NABA = processNABA()
names(NABA) = gsub("NABA_","",names(NABA))
NABA_SETS = names(NABA)

## Function to obtain ECM scores 
## seurat_obj = seurat object used to obtain scores 
## NABA_SETS = gene_sets used to obtain scores
## Output: matrix of scores
get_ECMscores = function(seurat_obj,NABA_SETS){
  
  DefaultAssay(seurat_obj) = "SCT"
  
  for(gset in NABA_SETS){
    features = list(gset = NABA[[gset]])
    #ctrl_genes = sum(NABA[[gset]] %in% rownames(seurat_obj))
    ctrl_genes = 35
    
    if(ctrl_genes>0){
      seurat_obj = AddModuleScore(object = seurat_obj, 
                                      features = features, 
                                      name = gset, 
                                      ctrl = ctrl_genes,
                                      seed = 77)
    }
  }
  
  mod_ids = paste(NABA_SETS,as.character(1),sep="")
  measured_sets = mod_ids[mod_ids %in% colnames(seurat_obj@meta.data)]

  ECMmat = seurat_obj@meta.data[,measured_sets]
  colnames(ECMmat) = gsub("1","",colnames(ECMmat))
  
  return(t(ECMmat))
  
}

# Immune cells data

ImmuneMarkers = read.csv(file = "./markers/ImmuneCellMarker.tsv",
                         sep = "\t", header = T,
                         stringsAsFactors = F)

ImmuneMarkers = ImmuneMarkers %>% dplyr::select("Cell.type","Gene") %>% 
  group_by(Cell.type) %>% summarise(tmp = list(Gene)) %>%
  deframe()

names(ImmuneMarkers) = gsub(pattern = "[ /-]",replacement = "_",names(ImmuneMarkers))
ImmuneSets = names(ImmuneMarkers)

get_ImmuneScores = function(seurat_obj,ImmuneSets,ImmuneMarkers){
  
  DefaultAssay(seurat_obj) = "SCT"
  ctrl_genes = min(unlist(lapply(ImmuneMarkers, length)))
  
  for(gset in ImmuneSets){
    features = list(gset = ImmuneMarkers[[gset]])

    if(ctrl_genes>0){
      seurat_obj = AddModuleScore(object = seurat_obj, 
                                  features = features, 
                                  name = gset, 
                                  ctrl = ctrl_genes,
                                  seed = 77)
    }
  }
  
  mod_ids = paste(ImmuneSets,as.character(1),sep="")
  measured_sets = mod_ids[mod_ids %in% colnames(seurat_obj@meta.data)]
  
  ECMmat = seurat_obj@meta.data[,measured_sets]
  colnames(ECMmat) = gsub("1","",colnames(ECMmat))
  
  return(t(ECMmat))
  
}