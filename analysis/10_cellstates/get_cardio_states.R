# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform a whole cell state characterization
#' 
#' 1) We identify genes that are representative for each state (40% cutoff)
#' 2) We check proportions of states in patients and keep states that are representative of a population (larger than 1% in 5 patients)
#' 3) After state filtering we perform linear mixed models to identify genes of interest
#' 4) We used those t-values for pathway activity estimation and tf activity estimation

library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(HDF5Array)
library(uwot)
library(lme4)
library(lmerTest)
library(progeny)
library(dorothea)
library(viper)
library(fgsea)
library(cowplot)
source("./analysis/utils/pseudobulk_utils.R")

# Load SCE object
group_variable <- "opt_state"
sc_data_file <- "./results/ct_data/cardiomyocyte/cardiomyocyte_states_sce/"
group_alias <- "state"
perc_thrsh <- 0.2
gene_filter_out <- "./results/ct_data/cardiomyocyte/gene_filter_df.rds"
n_samples_filt <- 5
state_prop_results <- "./results/ct_data/cardiomyocyte/state_proportion_kwtest.txt"
state_gene_list <- "./results/ct_data/cardiomyocyte/state_genelist.rds"
gsea_out <- "./results/ct_data/cardiomyocyte/state_gsea.txt"
progeny_out <- "./results/ct_data/cardiomyocyte/state_progeny.txt"
de_out <- "./results/ct_data/cardiomyocyte/state_de.txt"
tf_out <- "./results/ct_data/cardiomyocyte/state_tf.txt"
all_out <- "./results/ct_data/cardiomyocyte/states_funcomics.pdf"

sc_data <- loadHDF5SummarizedExperiment(sc_data_file)
sc_data[[group_variable]]  <- paste0(group_alias, sc_data[[group_variable]])
sc_meta <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")

sc_data_umap <- reducedDim(sc_data, "UMAP") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  left_join(sc_meta, by = "cell_id") 
  
# Here we need to get a matrix of 0 and 1 to know if a gene is expressed
# What is the minimum niche size?

n_cell_df <- colData(sc_data) %>% 
  as.data.frame() %>%
  group_by_at(group_variable) %>%
  summarise(n_spots = n())

min_spots <- n_cell_df %>%
  pull(n_spots) %>%
  min()

# First let's filter all genes that aree extremely lowly expressed
expressed_genes <- counts(sc_data) > 0
gene_ix <- rowSums(expressed_genes) > min_spots
sc_data <- sc_data[gene_ix, ]

# Then for each niche we will test specifically if the gene can be considered expressed
# Get cell_ids that belong to a class
niche_info <- colData(sc_data) %>% 
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  dplyr::select_at(c("cell_id", group_variable)) %>%
  group_by_at(group_variable) %>%
  nest() %>%
  rename("cell_ids" = data) %>%
  dplyr::mutate(cell_ids = map(cell_ids, ~ .x[[1]]))

niche_info <- niche_info %>%
  left_join(n_cell_df) %>%
  mutate(min_spots = (n_spots * perc_thrsh) %>% floor()) %>%
  dplyr::mutate(selected_genes = map2(cell_ids, min_spots, function(cids, mspots) {
    
    gene_ix <- rowSums(expressed_genes[, cids]) > mspots
    return(names(gene_ix[gene_ix]))
    
  }))

# This contains the filtered genes

#saveRDS(niche_info, gene_filter_out)

niche_info <- readRDS(gene_filter_out)

# Free some memory
rm(sc_data)
rm(expressed_genes)

# Part 2: Are all states representative of samples?
# This is to only consider states that aren't unique to a single
# sample

cell_state_prop <- sc_meta %>%
  dplyr::group_by_at(c("orig.ident", group_variable)) %>%
  summarise(n_cells_state = n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(n_cells_sample = sum(n_cells_state)) %>%
  ungroup() %>%
  mutate(props_state = n_cells_state/n_cells_sample)

state_counts <- cell_state_prop %>%
  dplyr::filter(props_state > 0.01) %>%
  group_by_at(group_variable) %>%
  summarize(n_samples = n())

filtered_states <- state_counts %>%
  dplyr::filter(n_samples >= n_samples_filt) %>%
  pull(group_variable)

cell_state_prop <- cell_state_prop %>%
  filter(get(group_variable) %in% filtered_states)

# N cells after filtering
n_cells_plot <- cell_state_prop %>%
  dplyr::select(orig.ident, n_cells_sample) %>%
  unique() %>%
  ggplot(aes(x = n_cells_sample, y = orig.ident)) +
  geom_bar(stat = "identity") +
  ylab("sample") +
  xlab("number of called cells")
  
# Prop of states per sample
sample_states_dots <- cell_state_prop %>%
  ggplot(aes(x = opt_state, y = orig.ident, size = props_state)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))

# UMAP
sc_data_umap_plt <- sc_data_umap %>%
  dplyr::filter(opt_state %in% filtered_states) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = opt_state)) +
  geom_point(size = 0.2) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size=4)))

# Part 3: Are these states changing between conditions?

sample_dict <- readRDS("./markers/snrna_patient_anns_revisions.rds") %>%
  mutate(batch = ifelse(grepl("CK3", sample_id), "B", "A"))

kwallis_res <- cell_state_prop %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id")) %>%
  group_by_at(group_variable) %>%
  nest() %>%
  mutate(kwallis_test = map(data, function(dat){
    broom::tidy(kruskal.test(props_state ~ patient_group, data = dat))
  })) %>%
  dplyr::select(kwallis_test) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value))

kwallis_res %>% 
    write.table(file = state_prop_results, quote = F, 
                col.names = T, row.names = F, sep = "\t")

kwallis_winners <- dplyr::filter(kwallis_res, p_corr < 0.15) %>%
  pull(opt_state)

# WATCH OUT MANUAL PLOT!!!
cell_state_prop_plt <- cell_state_prop %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id")) %>%
  mutate(opt_state = ifelse(opt_state %in% kwallis_winners, paste0(opt_state, "*"), opt_state)) %>%
  ggplot(aes(x = patient_group, y = props_state, color = patient_group)) +
  geom_boxplot() +
  facet_wrap(.~opt_state, nrow = 2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")

# Part 4: Pseudobulk analysis
pb_data <- readRDS("./results/ct_data/cardiomyocyte/ps_cardiomyocyte_states.rds")[[1]][["gex"]]
pb_meta <- colData(pb_data) %>% as.data.frame() %>% left_join(sample_dict, by = c("orig.ident" = "sample_id"))
pb_meta[, group_variable] <- paste0(group_alias, pb_meta[, group_variable])
pb_meta[, "col_id"] <- paste0(pb_meta$orig.ident, "_", pb_meta$opt_state)
pb_gex <- assay(pb_data)
colnames(pb_gex) <- pb_meta[, "col_id"]

# Filtering samples of interest
# Filtering niches from previous analysis and genes that are representative of a niche
niche_filter_ix <- pb_meta[, group_variable] %in% filtered_states
pb_gex <- pb_gex[niche_info$selected_genes %>% unlist() %>% unique(),
                 niche_filter_ix]

pb_meta <- pb_meta[niche_filter_ix, ]

# Then rejecting all pseudobulk profiles that come from less than one percent of population
pb_meta <- pb_meta %>%
  left_join(cell_state_prop %>%
              dplyr::select_at(c("orig.ident", group_variable, "props_state")),
            by = c("orig.ident", group_variable))

perc_filt <- which(pb_meta$props_state > 0.01)
pb_gex <- pb_gex[, perc_filt]
pb_meta <- pb_meta[perc_filt, ]

# Normalizing
pb_gex <- pb_gex %>%
  edgeR_filtering(expression_matrix = .,
                  min.count = 10, 
                  min.prop = 0.5) %>%
  cpm_norm(expression_matrix = .)

# UMAP - Why do I get patient effects?
# probably because we are not using the same genes as before

gex_umap <- umap(scale(t(pb_gex))) %>%
  as.data.frame()

colnames(gex_umap) <- c("UMAP1", "UMAP2")

rownames(gex_umap) <- pb_meta$col_id 

gex_umap %>%
  rownames_to_column("col_id") %>%
  left_join(pb_meta, by = "col_id") %>%
  ggplot(aes(x = UMAP1, y = UMAP2, 
             color = patient_group)) +
  geom_point(size = 2) +
  theme_classic()

gex_umap %>%
  rownames_to_column("col_id") %>%
  left_join(pb_meta, by = "col_id") %>%
  ggplot(aes(x = UMAP1, y = UMAP2, 
             color = opt_state)) +
  geom_point(size = 2) +
  theme_classic()

# Linear mixed models

pb_cpm_long <- pb_gex %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "col_id", values_to = "expr") %>%
  left_join(pb_meta, by = "col_id")

# pb_cpm_long <- pb_cpm_long %>%
#   group_by(gene) %>%
#   nest() %>%
#   mutate(mxeff_models = map(data, function(x) {
#     states <- x %>% pull(opt_state) %>% unique() %>% set_names()
#     state_column <- x %>% pull(opt_state)
#     
#     map(states, function(niche_id) { 
#       
#       niche_data <- x %>% 
#         dplyr::mutate(niche = ifelse(opt_state == niche_id, 
#                                      "test_niche",
#                                      "reference_niche")) %>%
#         dplyr::mutate(niche = factor(niche, 
#                                      levels = c("reference_niche",
#                                                 "test_niche")))
#       
#       fit <- lmerTest::lmer(expr ~ niche + (1 | orig.ident), data = niche_data)
#       
#     }) %>% enframe()
#     
#   })) 
# 
# test_resuls <- pb_cpm_long %>%
#   dplyr::select(gene, mxeff_models) %>%
#   unnest() %>%
#   mutate(model_coefficients = map(value, function(x) {
#     mem_res <- summary(x)
#     mem_res$coefficients  %>%
#       as.data.frame() %>% 
#       rownames_to_column("term") %>% 
#       pivot_longer(-term)
#   })) %>%
#   mutate(anova_pval = map(value, function(x) {
#     aov_res <- anova(x)
#     aov_res$`Pr(>F)`
#   }))
# 
# degs <- test_resuls %>%
#   dplyr::select(-value) %>%
#   unnest(anova_pval) %>%
#   dplyr::select(-model_coefficients)
# 
# degs_ext <- test_resuls %>%
#   dplyr::select(-c("value", "anova_pval")) %>%
#   dplyr::rename("state" = name) %>%
#   unnest(model_coefficients) %>%
#   dplyr::filter(term == "nichetest_niche",
#                 name %in% c("t value", "Pr(>|t|)")) %>%
#   dplyr::select(-term) %>%
#   mutate(name = ifelse(name == "t value", 
#                        "t_value", 
#                        "pval")) %>%
#   pivot_wider(names_from = name, 
#               values_from = value) %>%
#   ungroup() %>%
#   arrange(state,-t_value)
# 
# write.table(degs_ext, col.names = T, row.names = F, quote = F, sep = "\t",
#             file = "./results/ct_data/cardiomyocyte/cardio_state_mrks.txt")

degs_ext <- read.table("./results/ct_data/cardiomyocyte/cardio_state_mrks.txt",header = T,sep = "\t") %>%
  as_tibble()

plot_genes <- degs_ext %>% 
  group_by(state) %>%
  mutate(corr_pval = p.adjust(pval)) %>%
  dplyr::filter(corr_pval < 0.15) %>%
  slice(1:10) %>%
  pull(gene) %>%
  unique()

#Showing the expression of marker genes 
pb_mrkr_plt <- pb_cpm_long %>% 
  #dplyr::select(gene, data) %>%
  #unnest() %>%
  dplyr::select(gene, col_id, expr, opt_state) %>%
  dplyr::filter(gene %in% plot_genes) %>%
  arrange(gene, opt_state) %>%
  group_by(gene) %>%
  mutate(expr = (expr - mean(expr)) / sd(expr),
         gene = factor(gene,
                       levels = plot_genes),
         col_id = factor(col_id,
                         levels = col_id)) %>%
  ggplot(aes(x = gene, y = col_id, fill = expr)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2()

pb_mrkr_plt_ts <- degs_ext %>% 
  dplyr::filter(gene %in% plot_genes) %>%
  mutate(gene = factor(gene,
                       levels = plot_genes)) %>%
  ggplot(aes(y = gene, x = state, fill = t_value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2()

# Get representative genes for later:
# This is to enrich in visium slides

state_genes <- degs_ext %>% 
  arrange(state, -t_value) %>%
  dplyr::filter(pval <= 0.005, t_value > 0) %>%
  group_by(state) %>%
  dplyr::slice(1:200) %>%
  dplyr::select(gene, state) %>%
  nest() %>%
  deframe() %>%
  map(., ~.x[[1]])

saveRDS(state_genes, state_gene_list)

# Part 5. Funcomics

gene_sets <- readRDS(file = "./markers/Genesets_Dec19.rds")[["MSIGDB_CANONICAL"]]

# First GSEA
gsea_res <- degs_ext %>%
  dplyr::select(-pval) %>%
  group_by(state) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    set_names(x$t_value, x$gene) %>%
      na.omit()
  })) %>%
  rename("gls" = data) %>%
  mutate(gsea_stats = map(gls, function(t_vals) {
    
    fgsea(pathways = gene_sets,stats = t_vals)
  }))
  
gsea_res <- gsea_res %>%
  select(gsea_stats) %>%
  unnest() %>%
  dplyr::filter(padj < 0.15) %>%
  arrange(state,-(NES))

gsea_res %>% 
  dplyr::select(-leadingEdge) %>%
  write.table(., sep = "\t",
              row.names = F, col.names = T,
              quote = F, file = gsea_out)

# Then PROGENy
progeny_res <- degs_ext %>% dplyr::select(-pval) %>%
  pivot_wider(names_from = state,
              values_from = t_value) %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  progeny(.,scale = T,top = 500)

progeny_res %>%
  as.data.frame() %>%
  rownames_to_column("state") %>%
  pivot_longer(-state, names_to = "pathway", values_to = "progeny_scores") %>%
  write.table(., sep = "\t",
              row.names = F, col.names = T,
              quote = F, file = progeny_out)

progeny_plt <- progeny_res %>%
  as.data.frame() %>%
  rownames_to_column("state") %>%
  pivot_longer(-state, names_to = "pathway", values_to = "progeny_scores") %>%
  ggplot(aes(x = pathway, y = state, fill = progeny_scores)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Finally TF acts w/ viper and dorothea

min_targets <- 5

data(dorothea_hs, package = "dorothea")

regulons <- dorothea_hs %>% 
  dplyr::filter(confidence != "E")

regulons <- regulons %>% 
  dplyr::filter(target %in% (degs_ext$gene %>% unique()))

filtered_regulons <- regulons %>% 
  dplyr::select(tf, target) %>%
  group_by(tf) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  mutate(n_targets = map_dbl(data, length)) %>%
  dplyr::filter(n_targets >= 10) %>%
  dplyr::select(-n_targets) %>%
  deframe() %>%
  names()
  
regulons <- regulons %>%
  dplyr::filter(tf %in% filtered_regulons) %>%
  dorothea::df2regulon()

dorothea_res <- degs_ext %>%
  dplyr::select(-pval) %>%
  group_by(state) %>%
  nest() %>%
  mutate(data = map(data, function(x) {
    set_names(x$t_value, x$gene) %>%
      na.omit() %>%
      sort(decreasing = T)
  })) %>%
  rename("gls" = data) %>%
  mutate(viper_stats = map(gls, function(t_vals) {
    msviper(ges = t_vals, regulon = regulons, minsize  = 5)
  }))

dorothea_res <- dorothea_res %>%
  dplyr::select(viper_stats) %>%
  mutate(viper_stats = map(viper_stats, function(v_stats) {
    
    es_res <- v_stats$es
    TFs <- names(es_res$nes)
    tibble(TF = TFs,
           nes = es_res$nes[TFs],
           pvalue = es_res$p.value[TFs],
           size = es_res$size[TFs])
    
  })) %>%
  unnest() %>%
  mutate(padj = p.adjust(pvalue))

dorothea_res %>% arrange(state, -abs(nes)) %>%
  write.table(., sep = "\t",
              row.names = F, col.names = T,
              quote = F, file = tf_out)


dorothea_res_sign <- dorothea_res %>%
  arrange(state, -nes) %>%
  group_by(state) %>%
  slice(1:10) 

dorothea_plt <- dorothea_res %>%
  dplyr::select(state, TF, nes) %>%
  dplyr::filter(TF %in% dorothea_res_sign$TF) %>%
  dplyr::mutate(TF = factor(TF,
                            levels = unique(dorothea_res_sign$TF))) %>%
  ggplot(aes(y = TF, x = state, fill = nes)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  

# Final figure
# N cells, prop of cell-states, umap
upper_panel <- cowplot::plot_grid(n_cells_plot, sample_states_dots, 
                   sc_data_umap_plt, ncol = 3,
                   rel_widths = c(0.2, 0.4, 0.5), align = "hv")

lower_panel_left <- cowplot::plot_grid(cell_state_prop_plt, progeny_plt, ncol = 1, rel_heights = c(0.5,0.4))

lower_panel_right <- cowplot::plot_grid(pb_mrkr_plt_ts, dorothea_plt, nrow = 1)

lower_panel <- cowplot::plot_grid(lower_panel_left, lower_panel_right, nrow = 1, rel_widths = c(0.5, 0.4))


pdf(all_out, height = 18, width = 18)

plot_grid(upper_panel,lower_panel, ncol = 1, rel_heights = c(0.35,0.5))

dev.off()
