# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Put in order patient annotations

library(tidyverse)

# LOAD DATA ANNOTATIONS
# the objective is to map the visium_patients_ext annotation
# to the seeds that have a global_ID that allow us to merge them

visium_patients_ext <- read.table("./markers/visium_annotations_ext.txt", header = T, sep = "\t")
snrna_patients_seed <- read_csv("./markers/snrnaseq_patient_anns_revisions.csv")
visium_patients_seed <- read_csv("./markers/visium_patient_anns_revisions.csv")
new_anns <- read.table("./markers/visium_patient_anns_aryan.txt",sep = "\t", header = T) %>%
  dplyr::select(pipeline_names, novel.annotation.short) %>%
  dplyr::filter(pipeline_names != "") %>%
  dplyr::rename("region_novel" = novel.annotation.short)

# First, fix visium annotations
visium_patients_ext <- visium_patients_ext %>%
  mutate(visium = gsub(".+_", "", sample_id)) %>%
  select(sample_id, condition, region, patient_group, visium) %>%
  left_join(visium_patients_seed, by = "visium") %>%
  select(-visium) %>%
  left_join(new_anns, by = c("sample_id" = "pipeline_names"))

saveRDS(visium_patients_ext,
        file = "./markers/visium_patient_anns_revisions.rds")

# Second, transfer all the info, except sample_id to the snrnaseq

snrna_patients <- snrna_patients_seed %>% left_join(visium_patients_ext %>%
  select(-c("sample_id", "rep", "patient")) %>%
    unique(), by = "global_ID")

saveRDS(snrna_patients,
        file = "./markers/snrna_patient_anns_revisions.rds")

readRDS("./markers/snrna_patient_anns_revisions.rds") %>%
  write.table(col.names = T, row.names = F, quote = F, sep = ",", file = "./markers/snrna_patient_anns_revisions.csv")

# Third, create batch information for the rest of the samoles

visium_patients_ext <- readRDS("./markers/visium_patient_anns_revisions.rds")

visium_batch_info <- visium_patients_ext %>%
  dplyr::select(sample_id, patient) %>%
  dplyr::mutate(batch = ifelse(grepl("Visium", sample_id), "B", "A")) %>%
  dplyr::rename("orig.ident" = sample_id)  %>%
  write.table(col.names = T, 
              row.names = F, 
              quote = F, 
              sep = ",", 
              file = "./markers/visium_batch_ann.csv")

snrna_patients <- readRDS("./markers/snrna_patient_anns_revisions.rds")

snrna_batch_info <- snrna_patients %>%
  dplyr::select(sample_id, patient) %>%
  dplyr::mutate(batch = ifelse(grepl("CK3", sample_id), "B", "A")) %>%
  dplyr::rename("orig.ident" = sample_id)  %>%
  write.table(col.names = T, 
              row.names = F, 
              quote = F, 
              sep = ",", 
              file = "./markers/snrna_batch_ann.csv")
  




