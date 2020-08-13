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

## Function to group Dorothea regulons. 
## Input: A data frame containing Dorothea regulons, as stored in 
## https://github.com/saezlab/ConservedFootprints/tree/master/data
## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}

## Function to process SpaceRanger folders
## Input: Folder path
## Output: Seurat object ready for MISTy
process_visium = function(dir_path){
  ## <<Gene_expression>>
  dat = Load10X_Spatial(data.dir = dir_path)
  dat = SCTransform(dat, assay = "Spatial", verbose = FALSE)
  
  ## <<DOROTHEA>>
  dorothea_regulon_human = read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")
  # We obtain the regulons based on interactions with confidence level A, B and C
  regulon = dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%
    df2regulon()
  
  tf_act_mat = viper(eset = as.matrix(dat[["SCT"]]@data), 
                     regulon = regulon, nes = TRUE, 
                     method = "scale", minsize = 4, 
                     eset.filter = FALSE,
                     verbose = FALSE)
  
  ## Repeated regulon RFXAP
  tf_act_mat = tf_act_mat[!duplicated(tf_act_mat),]
  dat[['dorothea']] = CreateAssayObject(counts = tf_act_mat)
  
  
  ## <<PROGENY>>
  progeny_scores = progeny::progeny(expr = as.matrix(dat[["SCT"]]@data),
                                    scale=TRUE, 
                                    organism="Human", top=1000, perm=1)
  
  dat[['progeny']] = CreateAssayObject(counts = t(progeny_scores))

  ## <<Ligands>>
  
  opath_ligands = readRDS(file = "./markers/opath_ligands.rds")
  ligands = opath_ligands[opath_ligands$GeneSymbol %in% 
                          rownames(dat@assays$SCT@data),"GeneSymbol"] 
  
  ligands_gex = as.matrix(dat@assays$SCT@data[ligands,])
  ligands_expr = ligands_gex[rowSums(ligands_gex > 0)/ncol(ligands_gex) > .3,]
  
  dat[['ligands']] = CreateAssayObject(counts = ligands_expr)
  
  ## <<Heart markers scores>>
  markers_mat = readRDS("./markers/markers_matrix.rds")
  exprdata = as.matrix(dat[["SCT"]]@data)
  genes_in_markers = rownames(exprdata) %in% rownames(markers_mat)
  exprdata = exprdata[genes_in_markers,]
  markers_mat = markers_mat[rownames(exprdata),]
  
  ct_scores = scale(t(t(markers_mat) %*% exprdata))
  dat[['ctscores']] = CreateAssayObject(counts = t(ct_scores))
  
  ## <<ECM scores>>
  dat[['ECM']] = CreateAssayObject(counts = get_ECMscores(seurat_obj = dat,
                                                          NABA_SETS = NABA_SETS))
  
  return(dat)
}

#
#
#

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



