# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we estimate cell compositions from snRNAseq
library(tidyverse)

# Scell data meta --------------------------------------------------------------------------
scell_meta <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]][["annotations"]]
snrna_anns <- readRDS("./markers/snrna_patient_anns_revisions.rds")

# Generate cell type counts -----------------------------------------------------------
scell_props <- scell_meta %>%
  left_join(snrna_anns, by = c("orig.ident" = "sample_id")) %>%
  group_by(patient_id, cell_type) %>%
  summarize(sn_n_cells = length(patient_id)) %>%
  mutate(sn_prop_cells = sn_n_cells/sum(sn_n_cells)) %>%
  dplyr::mutate(cell_type = gsub("-","_", cell_type))

write.table(scell_props, file = "./results/compositions/snrna_compositions.txt", col.names = T, row.names = F, quote = F, sep = "\t")








