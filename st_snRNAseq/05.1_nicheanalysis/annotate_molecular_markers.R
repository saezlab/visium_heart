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

rm(all_niches)
#rm(meta_data)

# Generating markers per niche

colLabels(sc_data) <- factor(sc_data$molecular_niche)

wilcox_res <- scran::findMarkers(sc_data, test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrkrs <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest()

write_csv(all_mrkrs, file = "./processed_visium/integration/niche_mrkrs_mol_scran.csv")

# Defining useful niches
all_mrkrs <- read_csv("./processed_visium/integration/niche_mrkrs_mol_scran.csv")

sel_niches <- c("niche_0", 
                "niche_1", 
                "niche_3")

# Enrich functions in all markers
# MSigDB
gene_sets <- readRDS("./markers/Genesets_Dec19.rds")
hallmarks <- gene_sets$MSIGDB_HMARKS
canonical <- gene_sets$MSIGDB_CANONICAL

all_mrkrs_list <- all_mrkrs %>%
  dplyr::filter(summary.AUC > 0.55, FDR < 0.05) %>%
  dplyr::select(name, gene) %>%
  group_by(name) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

# First canonical pathways

ora_analysis_canonical <- map(all_mrkrs_list, GSE_analysis, Annotation_DB = canonical) %>%
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, p_value, corr_p_value, GenesInList)

write_csv(ora_analysis_canonical, "./results/niche_mapping/Spatial_snn_res.0.2/canonical_mrkr_ann.csv")

useful_gsets_canonical <- ora_analysis_canonical %>%
  dplyr::filter(corr_p_value < 0.15,
                name %in% sel_niches) %>%
  group_by(name) %>%
  dplyr::slice(1:5) %>%
  pull(gset) %>% unique()

ora_analysis_canonical_plts <- ora_analysis_canonical %>%
  mutate(logpval = -log10(corr_p_value)) %>%
  dplyr::filter(gset %in% useful_gsets_canonical,
                name %in% sel_niches) %>%
  dplyr::mutate(gset = factor(gset, levels = useful_gsets_canonical)) %>%
  ggplot(aes(y = gset,
             x = name,
             fill = logpval)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
  ylab("") +
  xlab("") +
  coord_equal()



pdf("./results/niche_mapping/Spatial_snn_res.0.2/canonical_mrkr_ann.pdf", height = 5.5, width = 15)

plot(ora_analysis_canonical_plts)

dev.off()

# Then hallmarks

ora_analysis_hmarks <- map(all_mrkrs_list, 
                           GSE_analysis, 
                           Annotation_DB = hallmarks) %>%
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, p_value, corr_p_value, GenesInList) 

write_csv(ora_analysis_hmarks, "./results/niche_mapping/Spatial_snn_res.0.2/hallmarks_mrkr_ann.csv")

useful_gsets_hmarks <- ora_analysis_hmarks %>%
  dplyr::filter(corr_p_value < 0.15,
                name %in% sel_niches) %>%
  group_by(name) %>%
  pull(gset) %>% unique()

ora_analysis_hmarks_plts <- ora_analysis_hmarks %>%
  mutate(logpval = -log10(corr_p_value)) %>%
  dplyr::filter(gset %in% useful_gsets_hmarks,
                name %in% sel_niches) %>%
  dplyr::mutate(gset = factor(gset, levels = useful_gsets_hmarks)) %>%
  ggplot(aes(y = gset,
             x = name,
             fill = logpval)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
  ylab("") +
  xlab("") +
  coord_equal()

pdf("./results/niche_mapping/Spatial_snn_res.0.2/hallmarks_mrkr_ann.pdf", height = 4, width = 11)

plot(ora_analysis_hmarks_plts)
write_csv()

dev.off()

# Show markers of useful niches
pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv") 

useful_genes <- all_mrkrs %>%
  dplyr::filter(summary.AUC > 0.55, FDR < 0.05,
                name %in% sel_niches) %>%
  group_by(name) %>%
  dplyr::slice(1:5) %>%
  pull(gene) %>%
  unique()




red_scdata <- logcounts(sc_data)[useful_genes %>% unique(), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  pivot_longer(-row_id,names_to = "gene",values_to = "log_gex")

red_scdata <- red_scdata %>% left_join(meta_data)

red_scdata <- left_join(red_scdata, pat_anns, by = c("orig.ident" = "sample_id"))

rm(sc_data)

gex_summary <- red_scdata %>%
  dplyr::filter(patient_group == "group_1",
                major_labl %in% c("BZ", "CTRL", "RZ"),
                molecular_niche %in% sel_niches) %>%
  group_by(molecular_niche, gene) %>%
  summarize(mean_gex = mean(log_gex)) %>%
  group_by(gene) %>%
  dplyr::mutate(std_mean_gex = (mean_gex - mean(mean_gex))/sd(mean_gex))

plt_genes <- gex_summary %>%
  dplyr::filter(gene %in% useful_genes,
                molecular_niche %in% sel_niches) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  ggplot(aes(x  = molecular_niche, y = gene, fill = std_mean_gex)) +
  geom_tile() +
  #geom_text() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()

plt_genes <- all_mrkrs %>%
  dplyr::filter(gene %in% useful_genes,
                name %in% sel_niches) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  mutate(sign = ifelse((FDR < 0.05 )& (summary.AUC > 0.55), "*","")) %>%
  ggplot(aes(x  = name, y = gene, fill = summary.AUC, label = sign)) +
  geom_tile() +
  geom_text() +
  scale_fill_gradient(low = "black", high = "yellow") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  coord_equal()
  
pdf("./results/niche_mapping/Spatial_snn_res.0.2/gene_mrkrs_ann.pdf", height = 4, width = 11)

plot(plt_genes)

write_csv(all_mrkrs %>%
            dplyr::filter(gene %in% useful_genes,
                          name %in% sel_niches), 
          "./results/niche_mapping/Spatial_snn_res.0.2/gene_mrkrs_ann.csv")

dev.off()




