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
patient_info <- readRDS("./markers/visium_patient_anns_revisions.rds")

# Raw qc stats
qc_feats <- c("sample_names",
              "Number of Spots Under Tissue")

slide_files <- paste0(path,
                      sample_names,
                      "/outs/metrics_summary.csv")

qc_stats <- tibble(sample_names,
                   qc_stats = map(slide_files, read_csv)) %>%
  unnest() %>%
  dplyr::select(all_of(qc_feats))

colnames(qc_stats) <- c("sample_names", "n_spots_under_tissue")

# Spot quantification
nuclei_quant_path <- "./results/Spotcounts/tissue_spot_counts_"

nuclei_data <- tibble(sample_name = sample_names,
                      nuclei_file = paste0(nuclei_quant_path, sample_names, ".csv")) %>%
  dplyr::mutate(nuclei_file = gsub("ium_", "ium",nuclei_file)) %>%
  dplyr::mutate(nuclei_quant = map(nuclei_file, read_csv)) %>%
  dplyr::select(-nuclei_file) %>%
  unnest() %>%
  dplyr::mutate(spot_id = paste0(sample_name, "..", barcode)) %>%
  dplyr::select(spot_id, count)

# Meta data from pseudobulk profiles
meta_data <- atlas_meta <- readRDS("./processed_visium/integration/ps_integrated_slides.rds")[[1]][["annotations"]] %>%
  rownames_to_column("spot_id") %>%
  dplyr::select(-patient) %>%
  dplyr::mutate(spot_id = strsplit(spot_id, "_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
  left_join(nuclei_data, by = "spot_id") %>%
  left_join(qc_stats, by = c("orig.ident" = "sample_names")) %>%
  left_join(patient_info, by = c("orig.ident" = "sample_id"))

# 

# Plot panels

ncells_p <- meta_data %>%
  group_by(patient_id) %>%
  summarise(ncells = (length(patient_id))) %>%
  ggplot(aes(x = patient_id, y = ncells)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of  \n spots")

nCount_RNA_p <- meta_data  %>%
  ggplot(aes(x = patient_id, y = log10(nCount_Spatial))) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of  \n UMIs (log10)")

nFeature_RNA_p <-  meta_data  %>%
  ggplot(aes(x = patient_id, y = log10(nFeature_Spatial))) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of  \n genes (log10)")

percent_mt_p <- meta_data  %>%
  ggplot(aes(x = patient_id, y = percent.mt)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Mitochondrial \n gene percentage")

nuclei_quant_p <- meta_data  %>%
  ggplot(aes(x = patient_id, y = count)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Quantified nuclei \n per spot")

# Generate panel
panel_p <- cowplot::plot_grid(ncells_p, nCount_RNA_p, 
                              nFeature_RNA_p, percent_mt_p,
                              nuclei_quant_p,
                              nrow = 5, ncol = 1, align = "hv")

pdf("./processed_visium/initial_qc/qc_distributions_spatial.pdf", height = 12, width = 6)
plot(panel_p)
dev.off()

# This is the table that we can use to correlate the loss of cells in snRNAseq

meta_data %>% 
  group_by(patient_id) %>%
  mutate(nspots = length(patient_id)) %>%
  dplyr::select(nspots, nCount_Spatial, nFeature_Spatial, count, percent.mt) %>%
  write_csv(.,
            file = "./processed_visium/initial_qc/qc_distributions_spatial.csv")

# This is the table that we can use to correlate the loss of cells in snRNAseq

# meta_data %>% 
#   group_by(patient_id) %>%
#   mutate(nspots = length(patient_id),
#          median_counts_spot = median(nCount_Spatial_filt),
#           median_genes_spot = median(nFeature_Spatial_filt)) %>%
#   dplyr::select(nspots, median_counts_spot, median_genes_spot, n_spots_under_tissue) %>%
#   unique() %>%
#   arrange(patient_id) %>%
#   mutate(n_spots_under_tissue = sum(n_spots_under_tissue)) %>%
#   unique() %>%
#   write.table(row.names = F, col.names = T, quote = F, sep = ",",
#             file = "./processed_visium/initial_qc/all_qcs_after_filtering.csv")
  


