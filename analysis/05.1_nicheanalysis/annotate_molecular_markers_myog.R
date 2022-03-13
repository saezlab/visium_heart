# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare the markers of niches coming
#' from similar structures

library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)
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

sel_niches <- c("niche_0", 
                "niche_1", 
                "niche_3")

sel_ids <- meta_data %>%
  dplyr::filter(molecular_niche %in% sel_niches) %>%
  pull(row_id)

sc_data <- sc_data[, sel_ids]

rm(all_niches)
#rm(meta_data)

# Generating markers per niche

colLabels(sc_data) <- factor(sc_data$molecular_niche)

wilcox_res <- scran::findMarkers(sc_data, test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrks_myog <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest()

write_csv(all_mrks_myog, file = "./processed_visium/integration/myogniche_mrkrs_mol_scran.csv")

# Comparison against all markers
all_mrks_myog <- read_csv("./processed_visium/integration/myogniche_mrkrs_mol_scran.csv")
all_mrkrs <- read_csv("./processed_visium/integration/niche_mrkrs_mol_scran.csv") %>%
  dplyr::rename("FDR_global" = FDR)

# All markers under the same threshold
# We keep the p-value because the population is
# the myogenic enriched group (upregulation within myogenic niches)
useful_all_mrks_myog <- all_mrks_myog %>%
  arrange(FDR) %>%
  dplyr::filter(summary.AUC > 0.5, FDR < 0.05,
                name %in% sel_niches) %>%
  arrange(name) %>%
  dplyr::select(name, gene, FDR, summary.AUC) %>%
  mutate(keep = T) %>%
  dplyr::rename("FDR_myogenic" = FDR,
                "summary.AUC_myogenic" = summary.AUC)

# I keep the AUC, because this ensures that
# this gene is representative of the niche

# I filter genes however that aren't the most expressed
# between myogenic niches (left join)

useful_all_mrkrs <- all_mrkrs %>%
  dplyr::select(-p.value) %>%
  arrange(FDR_global) %>%
  dplyr::filter(summary.AUC > 0.5, FDR_global < 0.05,
                name %in% sel_niches) %>%
  left_join(useful_all_mrks_myog, by = c("name", "gene")) %>%
  na.omit() %>%
  arrange(name, -summary.AUC)

useful_genes_v1 <- useful_all_mrkrs %>%
  group_by(name) %>%
  dplyr::slice(1:15) %>%
  pull(gene) %>%
  unique()

# Plot genes

useful_genes <- useful_genes_v1

red_scdata <- logcounts(sc_data)[useful_genes %>% unique(), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id,names_to = "gene",values_to = "log_gex")

red_scdata <- red_scdata %>% left_join(meta_data)

red_scdata <- left_join(red_scdata, pat_anns, by = c("orig.ident" = "sample_id"))

gex_summary <- red_scdata %>%
  dplyr::filter(patient_group == "group_1",
                major_labl %in% c("BZ", "CTRL", "RZ"),
                molecular_niche %in% sel_niches) %>%
  group_by(molecular_niche, gene) %>%
  summarize(mean_gex = mean(log_gex)) %>%
  group_by(gene) %>%
  dplyr::mutate(std_mean_gex = (mean_gex - mean(mean_gex))/sd(mean_gex)) 

gex_summary <- gex_summary %>%
  left_join(useful_all_mrkrs, by = c("molecular_niche" = "name", "gene")) %>%
  dplyr::mutate(sign = ifelse(FDR_myogenic < 0.05, "*", ""))

plt_genes <- gex_summary %>%
  dplyr::filter(gene %in% useful_genes,
                molecular_niche %in% sel_niches) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  ggplot(aes(x  = molecular_niche, y = gene, fill = std_mean_gex, label = sign)) +
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


pdf("./results/niche_mapping/Spatial_snn_res.0.2/gene_mrkrs_annv2.pdf", height = 10, width = 11)

plot(plt_genes)

write_csv(gex_summary,
          "./results/niche_mapping/Spatial_snn_res.0.2/gene_mrkrs_annv2.csv")

dev.off()

# Enrich CM marker genes
cm_state_mrkrs <- readRDS("./results/ct_data/state_genesets_list.rds")
cm_state_mrkrs <- cm_state_mrkrs[c("CM_healthy_CM", "CM_intermediate_CM", "CM_damaged_CM")] 

niche_marker_list <- useful_all_mrkrs %>%
  dplyr::select(name, gene) %>%
  group_by(name) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

map(niche_marker_list, GSE_analysis, Annotation_DB = cm_state_mrkrs)

map(cm_state_mrkrs, GSE_analysis, Annotation_DB = niche_marker_list)

# Plot marker genes

cm_state_mrkrs <- readRDS(file = "./results/ct_data/state_genesets.rds") %>%
  dplyr::filter(source %in% c("CM_healthy_CM", "CM_intermediate_CM", "CM_damaged_CM")) %>%
  dplyr::mutate(source = factor(source,
                                levels = c("CM_healthy_CM", "CM_intermediate_CM", "CM_damaged_CM"))) %>%
  arrange(source) %>%
  group_by(source) %>%
  dplyr::filter(target %in% useful_all_mrks_myog$gene) %>%
  dplyr::slice(1:5) %>%
  pull(target) %>%
  unique()

useful_genes <- cm_state_mrkrs

red_scdata <- logcounts(sc_data)[useful_genes %>% unique(), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id,names_to = "gene",values_to = "log_gex")

red_scdata <- red_scdata %>% left_join(meta_data)

red_scdata <- left_join(red_scdata, pat_anns, by = c("orig.ident" = "sample_id"))

gex_summary <- red_scdata %>%
  dplyr::filter(patient_group == "group_1",
                major_labl %in% c("BZ", "CTRL", "RZ"),
                molecular_niche %in% sel_niches) %>%
  group_by(molecular_niche, gene) %>%
  summarize(mean_gex = mean(log_gex)) %>%
  group_by(gene) %>%
  dplyr::mutate(std_mean_gex = (mean_gex - mean(mean_gex))/sd(mean_gex)) 

gex_summary <- gex_summary %>%
  left_join(useful_all_mrkrs, by = c("molecular_niche" = "name", "gene")) %>%
  dplyr::mutate(sign = ifelse(FDR_myogenic < 0.05, "*", ""))

plt_genes <- gex_summary %>%
  dplyr::filter(gene %in% useful_genes,
                molecular_niche %in% sel_niches) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  ggplot(aes(x  = molecular_niche, y = gene, fill = std_mean_gex, label = sign)) +
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

