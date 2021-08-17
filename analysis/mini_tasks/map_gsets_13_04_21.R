# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Marker expression
#' I am using the reduced states version of the atlas to reduce problems in interpretation
#' 
#' 

library(Seurat)
library(tidyverse)
library(stringr)
library(biomaRt)
source("./utils/funcomics.R")

#' Basic function to convert mouse genes to human
#' @param x: a vector of genes to be transformed 

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
convertmouse <- function(x){
  genes <- getLDS(attributes = c("mgi_symbol"), 
                  filters = "mgi_symbol", 
                  values = x, 
                  mart = mouse, 
                  attributesL = c("hgnc_symbol"),
                  martL = human, 
                  uniqueRows=T)
  
  return(genes)
}

# Main -------------------------------------------
atlas_dir <- "/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/integration/integrated_data_fordeconv.rds"
mi_atlas <- readRDS(atlas_dir)
gene_list <- list("MarkerAPM" = c("Nlrc5","Hsp90ab1","Ciita",
                                "Psmb10","Psmb9","Psme1",
                                "Psme2","B2m", "Tap1","Calr",
                                "Pdia3","Psme3","Rfx5", "Hsp90aa1", 
                                "Hsp90b1","Tapbpl", "Erap1"),
                  "MarkerAutophagy" = c("Tpp1","Cln3","Clta","Atp6v0c","Cd68","Naglu",
                                      "Hexb","Ctsa","Slc11a2","Smpd1","Tcirg1","Ap1s1",
                                      "Cd63","Galns","Lgmn","Gusb","Npc2","Atg7","Atg5",
                                      "Bag1","Npc1","Ap1m1", "Ap3d1","Gnptab","Mcoln1","Mfsd8",
                                      "Ids","Atp6v0b","Atg4d","Lamp2","Scarb2","Gib1"),
                  "MarkerSenescence" = c("Cdkn1a","Cdkn2a","App","Ctnnb1","Mapk1","Rac1","Arf1",
                                       "Junb","Bst1","E2f2","Il10","Il1b","Itgam",
                                       "Itgax","Lmnb1","Parp14","Tnf"))

gene_list <- map(gene_list, convertmouse)
gene_list <- map(gene_list, ~ .x[,"HGNC.symbol"])
#saveRDS(gene_list, file = "./markers/gene_list_12_04_21.rds")

# Get module score matrix for simplicity

mi_atlas <- getTF_matrix_MS(visium_slide = mi_atlas,
                                MS_regulon = gene_list,
                                assay = "RNA",
                                module_name = "gsets")


# Generalized violin ---------------------------------------

get_violins <- function(object, assay, group.by, gene_list) {
  
  vln_list <- map(gene_list, function(x) { 
    
    if(x %in% rownames(GetAssayData(object, 
                                    slot = "data", 
                                    assay = assay))) {
      
      
      plt <- Seurat::VlnPlot(object = object,
                      features = x,
                      assay = assay,
                      group.by = group.by)
      
      plt
      
    }
    
    })
  
}

# Get violins of all sets ----------------------------------------------------

ct_module_scores_vlns <- get_violins(object = mi_atlas,
                                  assay =  "gsets",
                                  group.by = "cell_type",
                                  gene_list = set_names(names(gene_list)))


pdf("./results/mini_tasks/violins_ms_cell_type.pdf")

walk(ct_module_scores_vlns, print)

dev.off()


cs_module_scores_vlns <- get_violins(object = mi_atlas,
                                     assay =  "gsets",
                                     group.by = "deconv_col",
                                     gene_list = set_names(names(gene_list)))

pdf("./results/mini_tasks/violins_ms_cell_states.pdf")

walk(ct_module_scores_vlns, print)

dev.off()

# Get violins of all genes ---------------------------------------------------

ct_genes_vlns <- map(gene_list, function(gset) {
  
  genes <- set_names(gset)
  
  get_violins(object = mi_atlas,
              assay =  "RNA",
              group.by = "cell_type",
              gene_list = genes)
})

walk(set_names(names(ct_genes_vlns)), function(x) {
  
  pdf(paste0("./results/mini_tasks/violins_cell_types_",x,".pdf"))
  
  walk(ct_genes_vlns[[x]], print)
  
  dev.off()
  
})

cs_genes_vlns <- map(gene_list, function(gset) {
  
  genes <- set_names(gset)
  
  get_violins(object = mi_atlas,
              assay =  "RNA",
              group.by = "deconv_col",
              gene_list = genes)
})


walk(set_names(names(ct_genes_vlns)), function(x) {
  
  pdf(paste0("./results/mini_tasks/violins_cell_states_",x,".pdf"))
  
  walk(cs_genes_vlns[[x]], print)
  
  dev.off()
  
})

















