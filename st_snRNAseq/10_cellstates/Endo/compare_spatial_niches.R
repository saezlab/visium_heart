# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare the markers of niches coming
#' from similar structures

library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)
library(Seurat)
library(ggpubr)
source("./analysis/utils/dea.R")

# visium data
sc_data <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")
sc_data <- scater::logNormCounts(sc_data)

# Patient annotation
pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

# Spot annotation
spot_anns <- read_csv("./results/niche_mapping/mol_clust_class.csv") %>%
  left_join(pat_anns, by = c("orig.ident" = "sample_id","patient_region_id"))

all_niches <- spot_anns %>%
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

colData(sc_data) <- meta_data %>%
  column_to_rownames("row_id") %>%
  DataFrame(.)

# Filter out spots not belonging to niche 9

spot_ix <- meta_data %>%
  dplyr::filter(molecular_niche == "niche_9")


# Check some QCs, like number of spots per patient

n_spots <- spot_ix %>%
  dplyr::group_by(patient_region_id) %>%
  summarise(n_spots = n()) %>%
  dplyr::arrange(-n_spots)

spot_ix <- spot_ix %>%
  pull(row_id)

# Filter data
sc_data <- sc_data[,spot_ix]

# Filter genes that are all 0
# Number of spots not containing 0
gene_sums = rowSums(counts(sc_data) != 0)

# Do you have at least 25% of the sampling spots?
gene_ix <- gene_sums >= (sum(n_spots$n_spots) * 0.25)

# Keep useful genes now

sc_data <- sc_data[gene_ix,]

# Now perform Wilcox using as identity the patient group

colLabels(sc_data) <- factor(sc_data$patient_group)

wilcox_res <- scran::findMarkers(sc_data, test.type = "wilcox", assay.type = "logcounts", pval.type = "all")

all_mrks <- lapply(wilcox_res, function(x) {
  
  x %>% 
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(-summary.AUC) %>%
    dplyr::select_at(c("gene", "p.value", "FDR", "summary.AUC"))
  
}) %>% enframe() %>%
  unnest()

selected_mrkrs <- all_mrks %>% dplyr::filter(summary.AUC > 0.55, FDR < 0.05)

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
  group_by(patient_group, gene) %>%
  summarize(mean_gex = mean(log_gex)) %>%
  group_by(gene) %>%
  dplyr::mutate(std_mean_gex = (mean_gex - mean(mean_gex))/sd(mean_gex))

gex_summary <- gex_summary %>%
  left_join(selected_mrkrs, by = c("patient_group" = "name", "gene")) %>%
  dplyr::mutate(sign = ifelse(FDR < 0.05, "*", ""))

plt_genes <- gex_summary %>%
  dplyr::filter(gene %in% useful_genes) %>%
  mutate(gene = factor(gene, levels = useful_genes)) %>%
  ggplot(aes(x  = patient_group, y = gene, fill = std_mean_gex, label = sign)) +
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

pdf("./results/cell_states/capillary_molniche9_degs.pdf")

plot(plt_genes)
write_csv(selected_mrkrs, "./results/cell_states/capillary_molniche9_degs_all.csv")

dev.off()

# Enrichment analysis of DEGs

mrkr_list <- selected_mrkrs %>%
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

pdf("./results/cell_states/capillary_molniche9_hmarks_gsa.pdf", width = 15, height = 5)

plot(plt_hmarks)
write_csv(gse_hmarks, "./results/cell_states/capillary_molniche9_hmarks_gsa.csv")

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

pdf("./results/cell_states/capillary_molniche9_canonical_gsa.pdf", width = 15, height = 5)

plot(plt_canonical)
write_csv(gse_canonical, "./results/cell_states/capillary_molniche9_canonical_gsa.csv")

dev.off()

# PROGENy

# # Get individual slide info ---------------------------------------------
# visium_folder = "./processed_visium/objects/"
# 
# visium_files <- list.files(visium_folder, full.names = F)
# visium_samples <- gsub("[.]rds", "", visium_files)
# 
# visium_df <- tibble(visium_file = paste0(visium_folder, 
#                                          visium_files),
#                     sample = visium_samples) %>%
#   mutate()
# 
# all_paths <- map(set_names(visium_df$visium_file, visium_df$sample), function(visium_file) { 
#   print(visium_file)
#   
#   path_score <- readRDS(visium_file) %>%
#     GetAssayData(., assay = "progeny") %>%
#     as.matrix() %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column("spot_id") %>%
#     pivot_longer(-spot_id)
#   
#   return(path_score)
# })
# 
# all_paths_red <- enframe(all_paths) %>%
#   dplyr::rename("orig.ident" = name) %>%
#   unnest() %>%
#   dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
#   left_join(spot_anns) %>%
#   dplyr::filter(spot_id %in% spot_ix)
# 
# 
# all_paths_red <- all_paths_red %>%
#   group_by(name, patient_group) %>%
#   dplyr::summarise(mean_act = mean(value)) %>%
#   dplyr::mutate(std_mean_act = (mean_act - mean(mean_act))/sd(mean_act))
# 
# 
# paths_mat <- all_paths_red %>%
#   dplyr::select(-mean_act) %>%
#   pivot_wider(names_from = name, values_from = std_mean_act) %>%
#   column_to_rownames("patient_group") %>%
#   as.matrix()
# 
# paths_plt <- ComplexHeatmap::Heatmap(t(paths_mat), name = "std_mean_act", cluster_columns = F)  
# 
# pdf("./results/cell_states/capillary_molniche9_progeny.pdf", width = 3.5, height = 5)
# 
# ComplexHeatmap::draw(paths_plt)
# 
# dev.off()
# 
# # Alternatively you can make boxplots
# my_comparisons <- list (c("group_1", "group_2"), c("group_1", "group_3"), c("group_2", "group_3"))
# 
# all_paths_red <- enframe(all_paths) %>%
#   dplyr::rename("orig.ident" = name) %>%
#   unnest() %>%
#   dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
#   left_join(spot_anns) %>%
#   dplyr::filter(spot_id %in% spot_ix)
# 
# all_paths_red <- all_paths_red %>%
#   group_by(name) %>%
#   nest() %>%
#   mutate(gplot = map2(name, data, function(nam, dat) {
#     
#     plt <- ggplot(dat,aes(x = patient_group, 
#                       color = patient_group,
#                       y = value)) +
#       geom_boxplot() +
#       theme_minimal() +
#       ggpubr::stat_compare_means(comparisons = my_comparisons) +
#       theme(axis.text.x = element_text(angle = 90,
#                                        hjust = 1,
#                                        vjust = 0.5),
#             axis.text = element_text(size = 10),
#             panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#       ylab("std.path. activity") +
#       ggtitle(nam)
#       
#     return(plt)
#     
#   }))
# 
# pdf("./results/cell_states/capillary_molniche9_progeny_bplots.pdf", height = 16, width = 14)
# 
# plot(cowplot::plot_grid(plotlist = all_paths_red$gplot, ncol = 4, align = "hv"))
# 
# dev.off()
# 
# # Summarize first by patient
# 
# all_paths_red <- enframe(all_paths) %>%
#   dplyr::rename("orig.ident" = name) %>%
#   unnest() %>%
#   dplyr::mutate(spot_id = paste0(orig.ident, "..", spot_id)) %>%
#   left_join(spot_anns) %>%
#   dplyr::filter(spot_id %in% spot_ix) %>%
#   group_by(name, patient_region_id, patient_group) %>%
#   summarize(mean_val = mean(value))
# 
# all_paths_red <- all_paths_red %>%
#   group_by(name) %>%
#   nest() %>%
#   mutate(pw_test = map(data, function(dat) {
#     
#     compare_means(mean_val ~ patient_group,  
#                   data = dat,
#                   method = "wilcox.test", 
#                   alternative = "two.sided") %>%
#       dplyr::select(group1, group2, p , p.adj)
#     
#   })) %>%
#   mutate(gplot = map2(pw_test, data, function(pw, dat) {
#     
#     max_val <- max(dat$mean_val) + 0.05
#     
#     plt <- ggplot(dat,aes(x = patient_group, 
#                           color = patient_group,
#                           y = mean_val)) +
#       geom_boxplot() +
#       geom_point() +
#       theme_minimal() +
#       ggpubr::stat_pvalue_manual(pw, label = "p.adj",
#                                  y.position = max_val, 
#                                  step.increase = 0.1,
#                                  tip.length = 0.01,size = 3) +
#       theme(axis.text.x = element_text(angle = 90,
#                                        hjust = 1,
#                                        vjust = 0.5),
#             axis.text = element_text(size = 10),
#             panel.border = element_rect(colour = "black", fill=NA, size=1)) +
#       ylab("std.path. activity") 
#     
#     return(plt)
#     
#   })) %>%
#   mutate(gplot = map2(name, gplot, function(nam,plt) {
#     
#     plt + ggtitle(nam)
#     
#   }))
# 
# pdf("./results/cell_states/capillary_molniche9_progeny_bplots_patient.pdf", height = 16, width = 14)
# 
# plot(cowplot::plot_grid(plotlist = all_paths_red$gplot, ncol = 4, align = "hv"))
# 
# dev.off()

# Running PROGENy per spot 

# We can do decoupleR
library(progeny)

network <- progeny::getModel(top = 1000) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "pathway", values_to = "likelihood") %>%
  dplyr::filter(likelihood != 0) %>%
  dplyr::mutate(mor = sign(likelihood)) %>%
  dplyr::mutate(likelihood =  abs(likelihood))


gset_mat <- decoupleR::run_wsum(network = network,
                                mat = logcounts(sc_data),
                                .source = "pathway",
                                .target = "gene",
                                .mor = "mor",
                                .likelihood = "likelihood",
                                times = 100)

write_csv(gset_mat, "./results/cell_states/capillary_molniche9_progeny_scores.csv")

progeny_res <- gset_mat %>%
  dplyr::filter(statistic == "norm_wsum") %>%
  arrange(-abs(score)) %>%
  left_join(meta_data %>% dplyr::select(row_id, patient_group, orig.ident ),
            by  = c("condition" = "row_id"))

progeny_res <- progeny_res %>%
  group_by(source, orig.ident, patient_group) %>%
  summarise(mean_score = mean(score))

my_comparisons <- list (c("group_1", "group_2"), c("group_1", "group_3"), c("group_2", "group_3"))

all_paths_red <- progeny_res %>%
  group_by(source) %>%
  nest() %>%
  mutate(pw_test = map(data, function(dat) {
    
    compare_means(mean_score ~ patient_group,  
                  data = dat,
                  method = "wilcox.test", 
                  alternative = "two.sided") %>%
      dplyr::select(group1, group2, p , p.adj)
    
  })) %>%
  mutate(gplot = map2(pw_test, data, function(pw, dat) {
    
    max_val <- max(dat$mean_score) + 0.05
    
    plt <- ggplot(dat,aes(x = patient_group, 
                          color = patient_group,
                          y = mean_score)) +
      geom_boxplot() +
      geom_point() +
      theme_minimal() +
      ggpubr::stat_pvalue_manual(pw, label = "p.adj",
                                 y.position = max_val, 
                                 step.increase = 0.1,
                                 tip.length = 0.01,size = 3) +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5),
            axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      ylab("std.path. activity") 
    
    return(plt)
    
  })) %>%
  mutate(gplot = map2(source, gplot, function(nam,plt) {
    
    plt + ggtitle(nam)
    
  }))


pdf("./results/cell_states/capillary_molniche9_progeny_bplots_patient.pdf", height = 16, width = 14)
 
 plot(cowplot::plot_grid(plotlist = all_paths_red$gplot, ncol = 4, align = "hv"))
 
dev.off()

