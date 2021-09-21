# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we will gfind if receptors have a gradient and in which samples?

library(tidyverse)
library(Seurat)

folder <- "./results/spark"

new_anns <- read.table("./markers/visium_patient_anns_aryan.txt",sep = "\t", header = T) %>%
  dplyr::select(pipeline_names, novel.annotation.short) %>%
  dplyr::filter(pipeline_names != "") %>%
  dplyr::rename("region_novel" = novel.annotation.short)

sample_dict <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  left_join(new_anns, by = c("sample_id" = "pipeline_names"))

ion_channels <- read_csv("./markers/ion_channels.csv")[[1]]

# Load spark results
spark_res <- tibble(files = list.files(folder, full.names = T)) %>%
  mutate(sample_id = map_chr(list.files(folder, full.names = F), 
                             ~ gsub("_spark.csv", "",.)),
         stat_res = map(files, ~ read_csv(.x))) %>%
  dplyr::select(sample_id, stat_res) %>%
  unnest() %>%
  dplyr::filter(gene %in% ion_channels)

# Here I filtered for ion receptors that are spatially variable
significant_ir <- spark_res %>% 
  dplyr::filter(adjustedPval < 0.15) %>%
  dplyr::filter() %>%
  dplyr::select(sample_id, gene) %>%
  left_join(sample_dict)

write.table(significant_ir, "./results/ion_channels/ionrec_sign_persample.txt",row.names = F, col.names = T, quote = F, sep = "\t")


# Which receptors have spatial patterns in several slides?
ion_rec_count <- significant_ir %>%
  group_by(gene) %>%
  summarize(n_slides = n()) %>%
  arrange(-n_slides)

write.table(ion_rec_count, "./results/ion_channels/ionrec_signcount.txt",row.names = F, col.names = T, quote = F, sep = "\t")

# Which slides have more genes - The ones from ischemic samples
sample_counts <- significant_ir %>%
  group_by(sample_id) %>%
  summarize(n_genes = n()) %>%
  arrange(-n_genes) %>%
  left_join(sample_dict)

pat_groups_test <- sample_counts$patient_group %>%
  unique() %>%
  set_names() %>%
  map(., function(x) {
    
    test_df <- sample_counts %>%
      dplyr::select(n_genes,patient_group) %>%
      mutate(patient_group = ifelse(patient_group == x, "test_group", "reference_group")) %>%
      mutate(patient_group = factor(patient_group,
                                    levels = c("test_group", "reference_group")))
    
    wilcox.test(n_genes ~ patient_group, test_df, alternative = "greater") %>%
      broom::tidy()
  }) %>%
  enframe() %>%
  unnest() %>%
  write.table("./results/ion_channels/wtest_nionrec.txt",row.names = F, col.names = T, quote = F, sep = "\t")

# Compare patient groups from visium
gene_per_sample <- significant_ir %>%
  group_by(sample_id) %>%
  nest() %>%
  mutate(data = map(data, ~.x[[1]]))

walk2(gene_per_sample$sample_id, gene_per_sample$data, function(sample_id, genes) {
  print(sample_id)
  
  file_path <- paste0("./processed_visium/objects/", sample_id, ".rds")
  out_path <- paste0("./results/ion_channels/plots/", sample_id, ".pdf")
  
  slide_visium <- readRDS(file_path)
  DefaultAssay(slide_visium) <- "SCT"
  
  pdf(out_path)
  
  walk(genes, function(g){
    plt <- SpatialFeaturePlot(slide_visium, features = g)
    plot(plt)
  })
  
  dev.off()
  
})
  











