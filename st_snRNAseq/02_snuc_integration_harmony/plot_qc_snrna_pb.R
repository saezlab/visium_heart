# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we rgenerate basic QC plots of an integrated object
#' These include:
#' -Number of cells
#' -Distribution of counts
#' -Distribution of features
#' -Distribution of percentage of mitochondrial genes

library(scater)
library(tidyverse)
library(cowplot)

# Folder with sample names
path <- "./visium_data/"
sample_names <- list.files(path)

# Patient info
patient_info <- readRDS("./markers/snrna_patient_anns_revisions.rds")

# Meta data from pseudobulk profiles
meta_data <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]][["annotations"]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(-patient) %>%
  dplyr::mutate(cell_id = strsplit(cell_id, "_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(cell_id = paste0(orig.ident, "..", cell_id)) %>%
  left_join(patient_info, by = c("orig.ident" = "sample_id"))

# Plot panels

ncells_p <- meta_data %>%
  group_by(patient_id) %>%
  summarise(ncells = (length(patient_id))) %>%
  ggplot(aes(x = patient_id, y = log10(ncells))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of \n cells (log10)")

nCount_RNA_p <- meta_data  %>%
  ggplot(aes(x = patient_id, y = log10(nCount_RNA))) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of  \n UMIs (log10)")

nFeature_RNA_p <-  meta_data  %>%
  ggplot(aes(x = patient_id, y = log10(nFeature_RNA))) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of  \n genes (log10)")


# Generate panel
panel_p <- cowplot::plot_grid(ncells_p, nCount_RNA_p, 
                              nFeature_RNA_p, 
                              nrow = 3, ncol = 1, align = "hv")

pdf("./processed_snrnaseq/initial_qc/qc_distributions_snrna.pdf", height = 8, width = 6)
plot(panel_p)
dev.off()

# Data set for reproducibility
meta_data %>%
  dplyr::select(patient_id,
                nFeature_RNA,
                nCount_RNA) %>%
  group_by(patient_id) %>%
  mutate(ncells = length(patient_id)) %>%
  write_csv(., file = "./processed_snrnaseq/initial_qc/qc_distributions_snrna.csv")

#meta_data %>% 
#  group_by(patient_id) %>%
#  summarise(ncells = length(patient_id),
#            median_counts_cell = median(nCount_RNA),
#            median_genes_cell = median(nFeature_RNA)) %>%
#  write.table(row.names = F, col.names = T, quote = F, sep = ",",
#              file = "./processed_snrnaseq/initial_qc/all_qcs_after_filtering_snrna.csv")


