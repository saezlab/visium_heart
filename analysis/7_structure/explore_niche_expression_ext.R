# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We explore in more detail the results from the mixed model
#' the objective is to find clearer descriptions of the niches

library(scater)
library(tidyverse)
library(philentropy)
library(corrplot)
library(edgeR)
library(ComplexHeatmap)
library(uwot)
library(lme4)
library(lmerTest)
library(fgsea)
library(compositions)
library(vsn)
source("./analysis/utils/pseudobulk_utils.R")

# Here we first load the original resuls ---------------------------------------

mixed_effects_res <- readRDS("./processed_visium/integration/mixed_effects_dea.rds")

mixed_effects_res$niche %>% unique() %>% saveRDS(., file = "./markers/main_niches.rds")

mixed_effects_res <- mixed_effects_res %>%
  group_by(niche) %>%
  dplyr::mutate(corr_pval = p.adjust(pval, "BH"))

# Generate object to get gene filtering ---------------------------------------

filtered_mixeff_res <- mixed_effects_res %>%
  dplyr::filter(corr_pval < 0.15)

# Checking how many unique genes ------------------------

gene_significance_counter <- filtered_mixeff_res$gene %>%
  table()

hist(gene_significance_counter)

# Getting the most unique genes

gene_cut <- names(gene_significance_counter[gene_significance_counter <= 30])

# Generate the dataframe that contains the marker genes
# p-value cut-off + exclusivity cutoff
niche_markers <- filtered_mixeff_res %>%
  dplyr::filter(gene %in% gene_cut) %>%
  dplyr::select(gene, niche) %>%
  nest("gene") %>%
  dplyr::rename(gene = data) %>%
  mutate(gene = map(gene, ~ .x[[1]])) %>%
  deframe()

saveRDS(niche_markers, "./markers/niche_markers_ps.rds")

# Filter t-values to include genes with significant difference and quite unique to the niche

mixed_effects_res <- mixed_effects_res %>%
  dplyr::filter(gene %in% gene_cut)
  
red_degs_ext <- mixed_effects_res %>% 
  dplyr::select(niche, gene, t_value) %>%
  pivot_wider(names_from = gene, values_from = t_value) %>%
  column_to_rownames("niche") %>%
  as.matrix()

draw(Heatmap(t((red_degs_ext)),
             show_row_names = FALSE,
             row_dend_side = "left",
             heatmap_legend_param = list(direction = "horizontal")),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")

# Let's plot the names of the genes of the top genes (10)

top_n <- 10

gene_cut <- filtered_mixeff_res %>%
  dplyr::filter(gene %in% gene_cut) %>%
  ungroup() %>%
  arrange(niche, -t_value) %>%
  group_by(niche) %>%
  dplyr::slice(1:top_n) %>%
  pull(gene)

red_degs_ext <- mixed_effects_res %>% 
  dplyr::select(niche, gene, t_value) %>%
  dplyr::filter(gene %in% gene_cut) %>%
  pivot_wider(names_from = gene, values_from = t_value) %>%
  column_to_rownames("niche") %>%
  as.matrix()


pdf("./results/niche_mapping/DE_genes_niches.pdf", height = 4.5, width = 16)
draw(Heatmap(red_degs_ext,
             heatmap_legend_param = list(direction = "horizontal")),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

# Enrich using ORA

# Function to do enrichment -----------------------------------------------------------------

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

# Main ----------------------------------------------------------------------------------

gene_sets <- readRDS(file = "./markers/Genesets_Dec19.rds")[["MSIGDB_CANONICAL"]]

niche_ora_results <- map(niche_markers, GSE_analysis, Annotation_DB = gene_sets)

niche_ora_results_filt <- niche_ora_results %>% 
  enframe() %>% 
  unnest() %>%
  dplyr::select(name, gset, corr_p_value) %>%
  dplyr::filter(corr_p_value < 0.05) %>%
  group_by(name) %>%
  slice(1:15)










