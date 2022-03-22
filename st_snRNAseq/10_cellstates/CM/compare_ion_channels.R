# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we check to which states
#' the ordered probes belong

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(HDF5Array)
library(viridis)
source("./analysis/utils/misty_pipeline.R")

CM_mrkrs <- readRDS("./cell_states/CM/annotation.rds")

# Generate the dictionary of channel sets ---------------------------------------------------

probes <- read_csv("./markers/ion_channels.csv",col_names = T) %>%
  dplyr::rename("channel" = `transmural ion channel genes`) %>%
  unique() %>%
  dplyr::pull(channel)

markers <- list("REACTOME_ION_CHANNEL_TRANSPORT" = readRDS("./markers/Genesets_Dec19.rds")$MSIGDB_CANONICAL[["REACTOME_ION_CHANNEL_TRANSPORT"]],
                "TRANSMURAL_ION_CHANNELS" = probes) %>%
  enframe(name = "gset",
          value = "gene") %>%
  unnest()

write_csv(markers, file = "./results/ion_plots/mrkr_analysis/ionchannel_gsets.csv")

CM_mrkrs %>%
  left_join(markers) %>%
  na.omit()

# Annotations ------------------------------------------------------

annotation_list <- readRDS("./processed_snrnaseq/cell_states/cellstate_annotation_list.rds")

all_annotations <- enframe(annotation_list) %>% 
  unnest() %>%
  dplyr::select(-name)

# Load scell data and annotate it
cms <- c("damaged_CM", "healthy_CM")
genes <- markers$gene %>% unique()

sc_data <- loadHDF5SummarizedExperiment("./results/ct_data/CM/CM_states_sce/")

meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("raw_id") %>%
  dplyr::mutate(spot_id =  strsplit(raw_id, "_") %>%
                  map_chr(., ~.x[[1]])) %>%
  left_join(all_annotations) %>%
  na.omit() %>%
  dplyr::filter(annotation %in% cms)

sc_data <- sc_data[, meta_data$raw_id]
sc_data$annotation <- meta_data$annotation
colLabels(sc_data) <- factor(sc_data$annotation)

# Keep only ion channel related genes
genes <- genes[genes %in% rownames(logcounts(sc_data))]

# Run Wilcoxon test ------------------------------------------------------

wilcox_res <- scran::findMarkers(sc_data[genes, ], test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrkrs <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest()

filt_mrkrs <- all_mrkrs %>%
  dplyr::filter(name == "healthy_CM") %>%
  dplyr::filter((summary.AUC >= 0.6) |(summary.AUC <= 0.4), FDR < 0.05) %>%
  left_join(markers) 

mrkr_order <- filt_mrkrs$gene


all_mrkrs %>%
  dplyr::filter(name == "damaged_CM") %>%
  dplyr::filter((summary.AUC >= 0.6) |(summary.AUC <= 0.4), FDR < 0.05) %>%
  left_join(markers) 

# Create barplot

fix.labels <- function(x){
  x + 0.5
}

brplot <- filt_mrkrs %>%
  dplyr::mutate(gene = factor(gene, levels = mrkr_order)) %>%
  ggplot(., aes(y = summary.AUC - 0.5, 
                x = gene,
                fill = gset)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = fix.labels) + ylab("AUC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  ggtitle("Differential expression of genes in healthy CMs")

pdf("./results/ion_plots/mrkr_analysis/top_mrkrs.pdf", height = 4, width = 5)

plot(brplot)

dev.off()


visium_slide <- readRDS("processed_visium/objects/AKK002_157781.rds") %>%
  positive_states(., assay = "cell_states") %>%
  filter_states(slide = .,
                by_prop = F,
                prop_thrsh = 0.10)

pdf("./results/ion_plots/mrkr_analysis/AKK002_15778_RYR2.pdf", height = 5, width = 5)

plot(SpatialFeaturePlot(visium_slide, features = "RYR2"))

dev.off()

pdf("./results/ion_plots/mrkr_analysis/AKK002_15778_ATP2B4.pdf", height = 5, width = 5)

plot(SpatialFeaturePlot(visium_slide, features = "ATP2B4"))

dev.off()

pdf("./results/ion_plots/mrkr_analysis/AKK002_15778_cmhealthy.pdf", height = 5, width = 5)

DefaultAssay(visium_slide) <- "cell_states_pos"
plot(SpatialFeaturePlot(visium_slide, features = "CM-healthy-CM") +
  scale_fill_viridis(option = "D"))

dev.off()


all_annotations %>%
  mutate(cell_id =  strsplit(raw_id, "_") %>%
           map_chr(., ~.x[[1]])) %>%
  write_csv(., file = "./processed_snrnaseq/cell_states/integrated_rnasamples_anns_wstates.csv")





