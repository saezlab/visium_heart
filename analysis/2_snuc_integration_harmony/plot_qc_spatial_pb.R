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
  dplyr::mutate(spot_id = strsplit(spot_id, "_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
  left_join(nuclei_data, by = "spot_id")

# Plot panels

ncells_p <- meta_data %>%
  group_by(orig.ident) %>%
  summarise(ncells = length(orig.ident)) %>%
  ggplot(aes(x = orig.ident, y = ncells)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Number of  \n spots")

nCount_RNA_p <- meta_data  %>%
  ggplot(aes(x = orig.ident, y = nCount_Spatial)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Number of  \n UMIs")

nFeature_RNA_p <-  meta_data  %>%
  ggplot(aes(x = orig.ident, y = nFeature_Spatial)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Number of  \n genes")

percent_mt_p <- meta_data  %>%
  ggplot(aes(x = orig.ident, y = percent.mt)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Mitochondrial \n gene percentage")

nuclei_quant_p <- meta_data  %>%
  ggplot(aes(x = orig.ident, y = count)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Quantified nuclei \n per spot")

# Generate panel
panel_p <- cowplot::plot_grid(ncells_p, nCount_RNA_p, 
                              nFeature_RNA_p, percent_mt_p,
                              nuclei_quant_p,
                              nrow = 5, ncol = 1, align = "hv")

pdf("./processed_visium/initial_qc/qc_distributions.pdf", height = 12, width = 6)
plot(panel_p)
dev.off()



meta_data %>% 
  group_by(orig.ident) %>%
  summarise(nspots = length(orig.ident),
            median_counts_spot = median(nCount_Spatial_filt),
            median_genes_spot = median(nFeature_Spatial_filt)) %>%
  print(n = 30) %>%
  write.table(row.names = F, col.names = T, quote = F, sep = ",",
            file = "./processed_visium/initial_qc/all_qcs_after_filtering.csv")
  


