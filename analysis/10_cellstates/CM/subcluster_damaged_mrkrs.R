# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we identify  markers in the integrated spatial data

library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)

sc_data <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")
sc_data <- scater::logNormCounts(sc_data)
cluster_info <- read_csv("./results/stressed_CMs/sup_figures/spot_ann.csv") %>%
  mutate(row_id = paste0(orig.ident, "..", spot_id)) %>%
  column_to_rownames("row_id")

# MSigDB
gene_sets <- readRDS("./markers/Genesets_Dec19.rds")
hallmarks <- gene_sets$MSIGDB_HMARKS
canonical <- gene_sets$MSIGDB_CANONICAL

# Generate mappable ids ---------------------------------------------------
meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  mutate(row_id = strsplit(row_id, "_") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(row_id = paste0(orig.ident, "..", row_id))

colnames(sc_data) <- meta_data$row_id

sc_data <- sc_data[, rownames(cluster_info)]

sc_data$cm_clust <- cluster_info[colnames(sc_data), "shared_label"]

colLabels(sc_data) <- factor(sc_data$cm_clust)

# Run Wilcoxon test ------------------------------------------------------

wilcox_res <- scran::findMarkers(sc_data, test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrkrs <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest() %>%
  dplyr::filter(summary.AUC > 0.55, FDR < 0.05)


all_mrkrs_list <- all_mrkrs %>%
  dplyr::select(name, gene) %>%
  group_by(name) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()


write_csv(all_mrkrs, file = "./results/stressed_CMs/sup_figures/clstr_mrkrs.csv")

# Enrich marker list

# Overrepresentation analysis

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
  
  ResultsDF = ResultsDF %>% 
    rownames_to_column("gset") %>% 
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"), 
              as.numeric) %>% 
    dplyr::arrange(corr_p_value,GenesInList)
  
  return(ResultsDF)
  
}

# Independent runs ----------------------------------

# Canonical pathways

ora_analysis_canonical <- map(all_mrkrs_list, GSE_analysis, Annotation_DB = canonical) %>%
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, p_value, corr_p_value, GenesInList)

write_csv(ora_analysis_canonical, "./results/stressed_CMs/sup_figures/canonical_chr.csv")

useful_gsets_canonical <- ora_analysis_canonical %>%
  dplyr::filter(corr_p_value < 0.15) %>%
     group_by(name) %>%
     dplyr::slice(1:5) %>%
  pull(gset) %>% unique()

ora_analysis_canonical_plts <- ora_analysis_canonical %>%
      mutate(logpval = -log10(corr_p_value)) %>%
  dplyr::filter(gset %in% useful_gsets_canonical) %>%
  dplyr::mutate(gset = factor(gset, levels = useful_gsets_canonical)) %>%
      ggplot(aes(y = gset,
                 x = name,
                 fill = logpval)) +
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
      ylab("") +
      xlab("") +
      coord_equal()

pdf("./results/stressed_CMs/sup_figures/canonical_chr.pdf", height = 4.5, width = 11)

plot(ora_analysis_canonical_plts)

dev.off()

# Hallmarks

ora_analysis_hmarks <- map(all_mrkrs_list, GSE_analysis, Annotation_DB = hallmarks) %>%
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 

write_csv(ora_analysis_hmarks, "./results/stressed_CMs/sup_figures/hmarks_chr.csv")

useful_gsets_hmarks <- ora_analysis_hmarks %>%
  dplyr::filter(corr_p_value < 0.15) %>%
  group_by(name) %>%
  dplyr::slice(1:5) %>%
  pull(gset) %>% unique()

ora_analysis_hmarks_plts <- ora_analysis_hmarks %>%
  mutate(logpval = -log10(corr_p_value)) %>%
  dplyr::filter(gset %in% useful_gsets_hmarks) %>%
  dplyr::mutate(gset = factor(gset, levels = useful_gsets_hmarks)) %>%
  ggplot(aes(y = gset,
             x = name,
             fill = logpval)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
  ylab("") +
  xlab("") +
  coord_equal()

pdf("./results/stressed_CMs/sup_figures/hmarks_chr.pdf", height = 4.5, width = 11)

plot(ora_analysis_hmarks_plts)

dev.off()

