# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Mapping of functional sets to SPARK results in single cells
#' 

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(ggplot2)

source("./visium_exploratory/slide_processing.R")

# Main

gene_sets = readRDS(file = "markers/Genesets_Dec19.rds")

spark_files = list.files("./results/spatial_cor/Archive")

slide_names = unlist(lapply(strsplit(spark_files,"_"), function(x) x[1]))
  
spark_files = set_names(paste0("./results/spatial_cor/Archive/",spark_files),
                        slide_names)

ora_res = lapply(spark_files, function(f){
  return(GSE_analysis(geneList = read.delim(file = f,sep = ",") %>%
                 dplyr::select(X,adjusted_pvalue) %>%
                 dplyr::filter(adjusted_pvalue <= 0.05) %>%
                 dplyr::select(X) %>% 
                 dplyr::filter(!grepl("^RP",X)) %>%
                   pull(),
               Annotation_DB = gene_sets$MSIGDB_CANONICAL))
}) %>% enframe("sample") %>% unnest() %>% 
  dplyr::filter(corr_p_value <= 0.05) %>%
  dplyr::arrange(sample,corr_p_value) %>%
  dplyr::group_by(sample) %>% unique()

spark_dots_df = ora_res %>%
  dplyr::slice(1:20) %>%
  ungroup() %>% group_by(gset) %>%
  mutate(ntimes = length(gset)) %>%
  ungroup() %>%
  dplyr::arrange(ntimes,sample) %>%
  dplyr::mutate(gset = substring(gset,1,50))
  
spark_dots = spark_dots_df %>%
  ggplot(aes(x=sample,y= factor(gset,
                                levels = unique(spark_dots_df$gset)),
             size = -log10(corr_p_value))) +
  geom_point() + theme_minimal() +
  theme(axis.text.y = element_text(size=8),
                       axis.text.x = element_text(size=8,angle = 90,
                                                  vjust = 0.5),
                       axis.title = element_blank(),
                       legend.text = element_text(size=8),
                       legend.title = element_text(size=8))


pdf(file = "results/spatial_cor/spark_dots_sl.pdf",
    width = 7, height = 7)

plot(spark_dots)

dev.off()

write.table(x = ora_res,
            file = "results/spatial_cor/spark_ora_sl.txt",
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")









