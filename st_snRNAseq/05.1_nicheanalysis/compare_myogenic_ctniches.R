# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare the markers of niches coming
#' from similar structures

library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)
library(decoupleR)
source("./analysis/utils/dea.R")

sc_data <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")
sc_data <- scater::logNormCounts(sc_data)
pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")


all_niches <- read_csv("./results/niche_mapping/niche_annotation_molct.csv") %>%
  dplyr::select_at(c("spot_id", "composition_niche", "Spatial_snn_res.0.2")) %>%
  dplyr::rename("molecular_niche" = `Spatial_snn_res.0.2`,
                "row_id" = spot_id) %>%
  mutate(composition_niche = paste0("niche_", composition_niche), 
         molecular_niche = paste0("niche_", molecular_niche))

meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  mutate(row_id = strsplit(row_id, "_") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(row_id = paste0(orig.ident, "..", row_id)) %>%
  left_join(all_niches)

colnames(sc_data) <- meta_data$row_id

sc_data$composition_niche <- meta_data$composition_niche
sc_data$molecular_niche <- meta_data$molecular_niche

# Just analyze data of the selected niches

sel_niches <- c("niche_8", 
                "niche_9", 
                "niche_7",
                "niche_1")

sel_ids <- meta_data %>%
  dplyr::filter(composition_niche %in% sel_niches) %>%
  pull(row_id)

sc_data <- sc_data[, sel_ids]

# Check some QCs, like number of spots per patient
spot_ix <-meta_data %>%
  dplyr::filter(composition_niche %in% sel_niches)

n_spots <- spot_ix %>%
  dplyr::group_by(composition_niche) %>%
  summarise(n_spots = n()) %>%
  dplyr::arrange(-n_spots)

# Filter genes that are all 0
# Number of spots not containing 0
gene_sums = rowSums(counts(sc_data) != 0)

# Do you have at least 25% of the sampling spots?
gene_ix <- gene_sums >= (sum(n_spots$n_spots) * 0.25)

# Keep useful genes now

sc_data <- sc_data[gene_ix,]


rm(all_niches)
#rm(meta_data)

# Generating markers per niche

colLabels(sc_data) <- factor(sc_data$composition_niche)

wilcox_res <- scran::findMarkers(sc_data, test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrks_myog <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest()

selected_mrkrs <- all_mrks_myog %>% dplyr::filter(summary.AUC > 0.55, FDR < 0.05)

write_csv(selected_mrkrs, "./results/niche_mapping/composition_niche/myogenicniches_degs_all.csv")

# Plotting differential expression

useful_genes <- selected_mrkrs %>%
  arrange(name, -summary.AUC) %>%
  group_by(name) %>%
  dplyr::slice(1:10) %>%
  pull(gene) %>%
  unique()

# Plot genes

red_scdata <- logcounts(sc_data)[useful_genes %>% unique(), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id,names_to = "gene",values_to = "log_gex")

red_scdata <- red_scdata %>% left_join(meta_data)

gex_summary <- red_scdata %>%
  group_by(composition_niche, gene) %>%
  summarize(mean_gex = mean(log_gex)) %>%
  group_by(gene) %>%
  dplyr::mutate(std_mean_gex = (mean_gex - mean(mean_gex))/sd(mean_gex))

gex_summary <- gex_summary %>%
  left_join(selected_mrkrs, by = c("composition_niche" = "name", "gene")) %>%
  dplyr::mutate(sign = ifelse(FDR < 0.05, "*", ""))

plt_genes <- gex_summary %>%
  dplyr::filter(gene %in% useful_genes) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  ggplot(aes(x  = composition_niche, y = gene, fill = std_mean_gex, label = sign)) +
  geom_tile() +
  geom_text(color = "darkgrey") +
  scale_fill_gradient2(low = ("blue"),
                       mid = "white",
                       high = ("red")) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()

pdf("./results/niche_mapping/composition_niche/myogenicniches_degs.pdf", width = 8, height = 3)

plot(plt_genes)

write_csv(gex_summary, "./results/niche_mapping/composition_niche/myogenicniches_degs.csv")


dev.off()

# Enrichment analysis of DEGs

mrkr_list <- all_mrks_myog %>% 
  dplyr::filter(summary.AUC > 0.5, FDR < 0.05) %>%
  group_by(name) %>%
  dplyr::select(gene) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

gsets <- readRDS("./markers/Genesets_Dec19.rds")

gse_hmarks <- map(mrkr_list, GSE_analysis , Annotation_DB = gsets$MSIGDB_HMARKS) %>%
  enframe() %>% unnest()

gse_canonical <- map(mrkr_list, GSE_analysis , Annotation_DB = gsets$MSIGDB_CANONICAL) %>%
  enframe() %>% unnest()

# filter and plot

prep_gsa_res <- function(gse_res) {
  
  useful_gsets <- gse_res %>%
    dplyr::filter(corr_p_value <= 0.05) %>%
    group_by(name) %>%
    dplyr::slice(1:10) %>%
    pull(gset) %>%
    unique()
  
  gse_red <- gse_res %>%
    dplyr::filter(gset %in% useful_gsets) %>%
    dplyr::mutate(gset = factor(gset, levels = useful_gsets),
                  logcorr_p_value = -log10(corr_p_value)) %>%
    dplyr::mutate(cap_p = ifelse(logcorr_p_value >=10, 10, logcorr_p_value))
  
  return(gse_red)
  
}

gse_hmarks_red <- prep_gsa_res(gse_res = gse_hmarks)

plt_hmarks <- ggplot(gse_hmarks_red, aes(x  = name, y = gset, fill = cap_p)) +
  geom_tile() +
  scale_fill_gradient(low = ("black"),
                      high = ("yellow")) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()

pdf("./results/niche_mapping/composition_niche/myogenicniches_hmarks_gsa.pdf", width = 15, height = 5)

plot(plt_hmarks)
write_csv(gse_hmarks, "./results/niche_mapping/composition_niche/myogenicniches_hmarks_gsa.csv")

dev.off()

gse_canonical_red <- prep_gsa_res(gse_res = gse_canonical)

plt_canonical <- ggplot(gse_canonical_red, aes(x  = name, y = gset, fill = cap_p)) +
  geom_tile() +
  scale_fill_gradient(low = ("black"),
                      high = ("yellow")) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()

pdf("./results/niche_mapping/composition_niche/myogenicniches_canonical_gsa.pdf", width = 15, height = 5)

plot(plt_canonical)
write_csv(gse_canonical, "./results/niche_mapping/composition_niche/myogenicniches_canonical_gsa.csv")
dev.off()

# We can do decoupleR of hypertrophic signature

network <- gsets$MSIGDB_CANONICAL %>%
  enframe() %>%
  unnest() %>%
  dplyr::filter(grepl("HYPERTROPHIC_CARDIOMYOPATHY",name)) %>%
  dplyr::mutate(mor = 1,
                likelihood = 1) %>%
dplyr::rename("target" = value, 
         "source" = name)

gset_mat <- decoupleR::run_wmean(network = network,
                      mat = logcounts(sc_data),
                      .source = "source",
                      .target = "target",
                      .mor = "mor",
                      .likelihood = "likelihood",
                      times = 100)

write_csv(gset_mat, "./results/niche_mapping/composition_niche/myogenicniches_HCM_wmean.csv")


norm_wmean_hcm <- gset_mat %>% 
  dplyr::filter(statistic == "norm_wmean") %>%
  dplyr::select(-c("statistic","p_value")) %>%
  left_join(spot_ix %>%
              dplyr::select(row_id, composition_niche),
            by = c("condition" = "row_id"))

#

my_comparisons <- list(c("niche_1", "niche_7"), c("niche_1", "niche_8"), c("niche_1", "niche_9"),
                        c("niche_7", "niche_8"), c("niche_7", "niche_9"), c("niche_8", "niche_9"))

hcm_boxplot <- ggplot(norm_wmean_hcm, aes(x = composition_niche, y = score)) +
  geom_violin() +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("normalized weighted mean") +
  ggpubr::stat_compare_means(comparisons = my_comparisons) +
  ggtitle("HYPER.CARDIOMYOPATHY")

pdf("./results/niche_mapping/composition_niche/myogenicniches_HCM_wmean.pdf", height = 4, width = 3.5)

plot(hcm_boxplot)

dev.off()








