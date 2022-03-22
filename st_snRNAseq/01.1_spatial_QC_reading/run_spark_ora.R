# Copyright (c) [2021] [Jovan Tanevski, Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we analyze the results of the SPARK analysis
#' folder
#' |
#' |---sample_spark.csv

library(tidyverse)
source("./analysis/utils/pca_utils.R")

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

folder <- "./results/spark"

sample_dict <- read.table("./markers/visium_annotations.txt",
                          sep = "\t", header = T)

sample_dict

# Load spark results

spark_res <- tibble(files = list.files(folder, full.names = T)) %>%
  mutate(sample_id = map_chr(list.files(folder, full.names = F), 
                             ~ gsub("_spark.csv", "",.)),
         stat_res = map(files, ~ read_csv(.x))) %>%
  dplyr::select(sample_id, stat_res) %>%
  unnest()

gene_sets <- readRDS(file = "./markers/Genesets_Dec19.rds")[["MSIGDB_CANONICAL"]]

# p-value filtering and ORA --------------------------------------------------------------------

spark_res <- spark_res %>%
  dplyr::filter(adjustedPval <= 0.001) %>%
  group_by(sample_id) %>%
  dplyr::select(gene) %>%
  nest("gene_list" = gene) %>%
  mutate(gene_list = map(gene_list, ~ .x[[1]])) %>%
  mutate(ora_results = map(gene_list, ~ GSE_analysis(geneList = .x,
                                                Annotation_DB = gene_sets))) %>%
  dplyr::select(ora_results) %>%
  unnest()

# writing results ------------------------------------------------------------------------------       

spark_res %>%
#dplyr::filter(corr_p_value < 0.001) %>%
  write_csv("./results/spark/spark_ora.csv", col_names = T)

# # let's use the p-values to cluster samples -----------------------------------------------------
# 
# # first get useful gsets
# 
# top_gsets <- spark_res %>% 
#   slice(1:30) %>%
#   ungroup() %>%
#   pull(gset) %>%
#   unique()
# 
# red_spark_res <- spark_res %>%
#   dplyr::filter(gset %in% top_gsets) %>%
#   ungroup() %>%
#   dplyr::select(gset, sample_id, corr_p_value) %>%
#   mutate(corr_p_value = -log10(corr_p_value)) %>%
#   pivot_wider(names_from = gset, values_from = corr_p_value) %>%
#   column_to_rownames("sample_id") %>%
#   as.matrix()
# 
# ora_pcs <- prcomp(x = red_spark_res[sample_dict$orig.ident,])
# 
# # print sample distribution
# 
# ggbiplot::ggbiplot(ora_pcs, obs.scale = 1, var.scale = 1, 
#                    labels = sample_dict$patient_sample, var.axes = F) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top')
# 
# # print relative importances
# 
# get_pc_importances(ora_pcs, pc = "PC2") %>% plot_pc_importances(.)
# 
# # Can we cluster results based on other sets (To be run if canonical pathways are used)
# 
# spark_res <- spark_res %>%
#   mutate(gset_short = substring(gset,1,50))
# 
# c_ann <- enframe(list(Meta = c("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
#                               "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY",
#                               "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHE",
#                               "REACTOME_COMPLEX_I_BIOGENESIS",
#                               "REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_C"),
#                      
#                      Cardiomyocyte_function = c("KEGG_CARDIAC_MUSCLE_CONTRACTION",
#                                                 "REACTOME_MUSCLE_CONTRACTION",
#                                                 "KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM",
#                                                 "KEGG_DILATED_CARDIOMYOPATHY"),
#                      
#                      ECM = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
#                              "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
#                              "NABA_CORE_MATRISOME",
#                              "NABA_ECM_GLYCOPROTEINS",
#                              "REACTOME_ECM_PROTEOGLYCANS",
#                              "NABA_MATRISOME"),
#                      
#                      Immune_process = c("REACTOME_NEUTROPHIL_DEGRANULATION",
#                                         "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
#                                         "REACTOME_INNATE_IMMUNE_SYSTEM"),
#                      
#                      Signaling = c("REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
#                                    "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES",
#                                    "REACTOME_PROGRAMMED_CELL_DEATH",
#                                    "PID_PDGFRB_PATHWAY",
#                                    "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
#                                    "REACTOME_SIGNALING_BY_INTERLEUKINS",
#                                    "REACTOME_VESICLE_MEDIATED_TRANSPORT",
#                                    "REACTOME_RHO_GTPASES_ACTIVATE_ROCKS")),name = "annotation",value = "gset_short") %>%
#   unnest() %>% 
#   left_join(data.frame("annotation" = c("Meta","Cardiomyocyte_function","ECM","Immune_process","Signaling"),
#                        "color" = c("darkblue","darkgreen","darkred","purple","darkorange")))
# 
# 
# # Creating plot
# 
# hmap_df = spark_res %>%
#   dplyr::select(gset_short, sample_id, corr_p_value) %>%
#   arrange(-corr_p_value) %>%
#   mutate(gset_short = factor(gset_short,levels = c_ann$gset_short),
#          sample_id = factor(sample_id, levels = c("157771", "157782", "157785", 
#                                             "157772", "157777", "157781", 
#                                             "157779", "157775")))
# 
# hmapsparkv2 = hmap_df %>%
#   dplyr::mutate(vis_p_value = -log10(corr_p_value)) %>%
#   dplyr::mutate(vis_p_value = ifelse(vis_p_value >=20,
#                                      20, vis_p_value)) %>%
#   ggplot(aes(x = sample_id, 
#              y = gset_short, 
#              fill = vis_p_value)) + 
#   geom_tile(na.rm = T) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90,vjust=1),
#         axis.text.y = element_text(hjust = 1, colour = c_ann$color)) +
#   scale_fill_gradient(
#     low = "black",
#     high = "yellow")
# 
# 
