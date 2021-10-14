# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will add cell-death signatures to the slides
#' 
#' Our assumption is that in all the IZ slides we will have a greater expression

library(Seurat)
library(tidyverse)
source("./analysis/utils/funcomics.R")

# First, prove that slides with more loss of spots
# have less cells

snrna_qc <- read_csv("./processed_snrnaseq/initial_qc/all_qcs_after_filtering_snrna.csv")
spatial_qc <- read_csv("./processed_visium/initial_qc/all_qcs_after_filtering.csv")

pat_meta <- left_join(snrna_qc, spatial_qc) %>%
  dplyr::mutate(perc_spot_loss = 1 - (nspots/n_spots_under_tissue))

patient_ann  <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  left_join(pat_meta)

pdf("results/cell_death/spot_loss_ncells.pdf", height = 4, width = 4)
print(ggplot(pat_meta, aes(x = log10(ncells), y = perc_spot_loss)) + geom_point())
dev.off()

# First find gene sets of interest
gsets <- readRDS("./markers/Genesets_Dec19.rds")
gsets <- c(gsets$MSIGDB_HMARKS, gsets$MSIGDB_CANONICAL)
  
apoptosis <- gsets[grepl("apoptosis", 
                         names(gsets),
                         ignore.case = T)]
necrosis <- gsets[grepl("necrosis", 
                        names(gsets),
                        ignore.case = T)]

death_gsets <- c(apoptosis, necrosis)

# Add these genesets into filtered slides...

visium_folder <- "./processed_visium/objects/"
sample_file <- list.files(visium_folder)

param_df <- tibble(sample_id = gsub("[.]rds", "", 
                                    sample_file),
                   slide_file = paste0(visium_folder,
                                       sample_file))


estimate_death_scores <- function(slide_file) {
  
  visium_slide <- readRDS(slide_file)
  
  visium_slide <- getTF_matrix_MS(visium_slide = visium_slide,
                                  MS_regulon = death_gsets,
                                  assay = "SCT",
                                  module_name = "death_gsets")
  
  sum_vector <- GetAssayData(visium_slide, assay = "death_gsets") %>%
    rowSums() 
  
  mean_vector <- GetAssayData(visium_slide, assay = "death_gsets") %>%
    rowMeans() 
  
  tibble("gset" = names(sum_vector), 
         "sum_score" = sum_vector,
         "mean_score" = mean_vector)
  
}

cell_death_scores <- param_df %>%
  mutate(death_scores = map(slide_file, estimate_death_scores)) %>%
  dplyr::select(sample_id, death_scores) %>%
  unnest()

selected_gsets <- c("HALLMARK-APOPTOSIS")

plot_data <- cell_death_scores %>%
  left_join(patient_ann) %>%
  dplyr::filter(gset %in% selected_gsets) %>%
  group_by(gset, patient_id) %>%
  mutate(sum_score = mean(sum_score),
         mean_score = mean(mean_score)) %>%
  dplyr::select(sum_score, mean_score, ncells, perc_spot_loss, patient_group) %>%
  dplyr::mutate(basic_def = ifelse(patient_group == "group_1", "high", "low")) %>%
  unique()

wilcox.test(ncells ~ basic_def, plot_data, alternative = "greater")
wilcox.test(mean_score ~ basic_def, plot_data, alternative = "less")

ncells_plot <- plot_data %>%
  ggplot(aes(x = patient_group, y = ncells)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size=12),
        axis.text.x = element_text(angle = 90)) +
  ylab("N nuclei in snRNAseq") +
  xlab("")


apoptosis_plot <- plot_data %>%
  ggplot(aes(x = patient_group, y = mean_score)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size=12),
        axis.text.x = element_text(angle = 90)) +
  ylab("Apoptosis signature") +
  xlab("")


pdf("results/cell_death/cell_death_qc.pdf", height = 3.5, width = 5)

cell_death_qc <- cowplot::plot_grid(ncells_plot, apoptosis_plot, align = "hv")
plot(cell_death_qc)

dev.off()




