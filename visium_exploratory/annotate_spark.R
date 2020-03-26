# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Mapping of functional sets to SPARK results

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)

## Function to perform Hypergeometric Tests for gene set enrichment 
## Input:
## geneList = query gene list to enrich
## Annotation_DB = a list of gene sets to enrich in geneList
## Output: data frame with results
GSE_analysis = function(geneList,Annotation_DB){
  
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
  
  return(ResultsDF)
  
}

#' Recover condition specific responsive genes
#' from spatial autocorrelation results
#' It always takes compares to Common results
get_condgenes = function(cor_res){
  
  # Creating iterator of conditions
  cnd_groups = unique(cor_res$Group)
  names(cnd_groups) = cnd_groups
  
  cnd_genes = lapply(cnd_groups, function(x){
    cor_res %>% filter(Group == x) %>% select(Gene) %>% pull()
  })
  
  
  cnd_genes[c("Borderzone","Chronic",
              "Healthy","MI")] = lapply(cnd_genes[c("Borderzone","Chronic",
                                                    "Healthy","MI")], function(x){
                                                      setdiff(x, cnd_genes$Common)
                                                    })
  return(cnd_genes)
  
}

#' Recover condition specific responsive genes
#' from spatial autocorrelation results
#' It always takes compares to Common results
annotate_autocor = function(cnd_genes,out_name){
  
  canonical_results = lapply(cnd_genes, GSE_analysis, 
                             Annotation_DB = gene_sets$MSIGDB_CANONICAL)
  canonical_results = lapply(canonical_results, function(x) x %>% 
                               rownames_to_column("gset") %>% 
                               filter(as.numeric(corr_p_value)<0.1))  %>% 
    enframe("Condition_group") %>% 
    unnest() %>% mutate_at(c("GenesInPathway","GenesInList",
                             "p_value","corr_p_value"), 
                           as.numeric) %>% 
    dplyr::arrange(Condition_group,corr_p_value)
  
  write.table(canonical_results,quote = F,sep = "\t",col.names = T,
              row.names = F, file = sprintf("results/%s_annotation_canonical.txt",
                                            out_name))
  
  ct_results = lapply(cnd_genes, GSE_analysis, Annotation_DB = cts_genes)
  ct_results = lapply(ct_results, function(x) x %>% 
                        rownames_to_column("gset") %>% 
                        filter(as.numeric(corr_p_value)<0.1)) %>% enframe("Condition_group") %>% 
    unnest() %>% mutate_at(c("GenesInPathway","GenesInList",
                             "p_value","corr_p_value"), 
                           as.numeric) %>% dplyr::arrange(Condition_group,corr_p_value)
  
  write.table(ct_results,quote = F,sep = "\t",col.names = T,
              row.names = F, file = sprintf("results/%s_annotation_celltypes.txt",
                                            out_name))
  
  return(list("canonical" = canonical_results,
         "cell_type" = ct_results))
  
}

# MAIN #
gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")
spark_res = read.csv("./results/spark_intersections.csv", header = T,
                     stringsAsFactors = F)
moran_res = read.csv("./results/moran_intersections.csv", header = T,
                     stringsAsFactors = F)
ct_markers = readRDS(file = "markers/wilcox_res.rds") %>%
  dplyr::group_by(cluster) %>% dplyr::arrange(p_val) %>%
  dplyr::slice(1:200) %>% select(gene) %>% ungroup()


# Creating iterator
cts = unique(as.character(ct_markers$cluster))
names(cts) = cts
# Getting cell type markers list
cts_genes = lapply(cts, function(x){
  ct_markers %>% filter(cluster == x) %>% select(gene) %>% pull()
})

# Annotate

spark_ann = annotate_autocor(cnd_genes = get_condgenes(spark_res),
                 out_name = "spark")

moran_ann = annotate_autocor(cnd_genes = get_condgenes(moran_res),
                             out_name = "moran")











