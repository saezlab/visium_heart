# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script I find myeloid-fibro probes for RNAscope
#' 
#' Fibro markers: Markers of state 6 and 5 that are not overlapping, highly expressed in macrophages
#' Myeloid markers: Markers of 0 and 5 that are not overlapping, highly expressed in macrophages

library(tidyverse)
library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(scuttle)
library(Seurat)
source("./analysis/utils/pseudobulk_utils.R")

# Patient_annotations

pat_ann <- readRDS("./markers/visium_patient_anns_revisions.rds")
pat_ann <- pat_ann %>%
  dplyr::filter(patient_group == "group_2")

# Load compositions of all slides here:
c2l_comp <- readRDS("./markers/comp2cl.rds")

# Load all slide info
int_slides <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")

meta_data_slides <- colData(int_slides) %>%
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  dplyr::mutate(spot_id = strsplit(spot_id, "_") %>%
                  map_chr(., ~.x[[1]]))

# Filter to only have ischemic patients
p_ix <- which(meta_data_slides$orig.ident %in% pat_ann$sample_id)
int_slides <- int_slides[, p_ix]
meta_data_slides <- meta_data_slides[p_ix,]

# Create a mask of myeloid = 0.175, fibro = 0.175
spots_interest <- c2l_comp %>%
  dplyr::filter(name %in% c("Myeloid")) %>%
  dplyr::mutate(cell_present = ifelse(value >= 0.175,
                                      1, 0)) %>%
  group_by(sample_id, spot_id) %>%
  summarise(ncells = sum(cell_present)) %>%
  ungroup()

# Add this info to meta_data and filter spots that aren't useful
meta_data_slides <- left_join(meta_data_slides, spots_interest,
                              by = c("orig.ident" = "sample_id",
                                     "spot_id"))

f_ix <- which(meta_data_slides$ncells > 0)
int_slides <- int_slides[, f_ix]
meta_data_slides <- meta_data_slides[f_ix,]

# How many spots per slide
nspots <- meta_data_slides %>% 
  group_by(orig.ident) %>% 
  summarise(n_use_spots = n())

# Filter info from slides with less than a 100 spots

useful_slides <- nspots %>%
  dplyr::filter(n_use_spots >= 100) %>%
  pull(orig.ident )

s_ix <- which(meta_data_slides$orig.ident %in% useful_slides)
int_slides <- int_slides[, s_ix]
meta_data_slides <- meta_data_slides[s_ix,]

# Identify exclusive marker genes of fibroblast state
fib_selected_states <- c("state5", "state6")
fib_markers <- read_table2("./results/ct_data/Fib/state_mrks.txt") %>%
  dplyr::filter(state %in% fib_selected_states,
                FDR < 0.15,
                logFC > 0) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  group_by(state) %>%
  nest()

state5_genes <- fib_markers %>%
  dplyr::filter(state == "state5") %>%
  pull(data) %>%
  enframe() %>%
  unnest() %>%
  pull(gene)
  
state6_genes <- fib_markers %>%
  dplyr::filter(state == "state6") %>%
  pull(data) %>%
  enframe() %>%
  unnest() %>%
  pull(gene)

f_state5_genes <- base::setdiff(state5_genes, state6_genes) %>%
  enframe(value = "gene") %>%
  mutate(name = "f_state5")

f_state6_genes <- base::setdiff(state6_genes, state5_genes) %>%
  enframe(value = "gene") %>%
  mutate(name = "f_state6")

# Identify exclusive marker genes of Myeloid state
myeloid_selected_states <- c("state0", "state5")
myeloid_markers <- read_table2("./results/ct_data/Myeloid/state_mrks.txt") %>%
  dplyr::filter(state %in% myeloid_selected_states,
                FDR < 0.15,
                logFC > 0) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  group_by(state) %>%
  nest()
  
m_state5_genes <- myeloid_markers %>%
  dplyr::filter(state == "state5") %>%
  pull(data) %>%
  enframe() %>%
  unnest() %>%
  pull(gene)

m_state0_genes <- myeloid_markers %>%
  dplyr::filter(state == "state0") %>%
  pull(data) %>%
  enframe() %>%
  unnest() %>%
  pull(gene)

m_state5_genes <- base::setdiff(m_state5_genes, m_state0_genes) %>%
  enframe(value = "gene") %>%
  mutate(name = "m_state5")

m_state0_genes <- base::setdiff(m_state0_genes, m_state5_genes) %>%
  enframe(value = "gene") %>%
  mutate(name = "m_state0")

# Ok, put all genes together

all_candidates <- bind_rows(f_state5_genes, f_state6_genes,
                            m_state5_genes, m_state0_genes)

all_candidates <- all_candidates %>%
  dplyr::filter(gene %in% rownames(int_slides))

candidate_genes <- unique(all_candidates$gene)

# Now the interesting part
int_slides <- int_slides[candidate_genes, ]
int_slides@assays@data$logcounts <- cpm_norm(expression_matrix = counts(int_slides))

gene_expression <- t(logcounts(int_slides)) > 0 %>%
  as.matrix()

gene_expreession_summ <- colSums(gene_expression)

gene_expreession_summ <- sort(gene_expreession_summ, decreasing = T)

gene_expression_df <- gene_expression %>%
  as.data.frame()
  
expression_df <- bind_cols(meta_data_slides[, c("spot_id", "orig.ident")], 
                           gene_expression_df) %>%
  pivot_longer(-c("spot_id", "orig.ident"))

# Now try to get genes with high expression in every slide

prop_expr_by_slide <- expression_df %>%
  group_by(orig.ident) %>%
  nest() %>%
  mutate(data = map(data, function(dat) {
    
    dat %>%
      group_by(name) %>%
      summarize(n_spots = n(),
             n_expr = sum(value)) %>%
      dplyr::mutate(prop_expr = n_expr/n_spots)
    
  })) %>%
  unnest() %>%
  arrange(name)

mean_prop <- prop_expr_by_slide %>%
  ungroup() %>%
  group_by(name) %>%
  summarise(mean_prop = mean(prop_expr)) %>%
  arrange(- mean_prop)

# Let's start with a conservative cutoff of 70%
genes_to_cor <- mean_prop %>%
  dplyr::filter(mean_prop >= 0.4) %>%
  pull(name)

# Now that we have our genes to correlate and our spots to correlate
rm(c2l_comp)
rm(expression_df)
rm(gene_expression)
rm(gene_expression_df)
rm(int_slides)
rm(log_counts)

# Get correlations between states and genes

get_gene_state_cor <- function(dat, f_path) {
  
  visium_slide <- readRDS(f_path)[, dat$spot_id]
  
  sel_states <- c("gene", "Fib-state5", "Fib-state6",
                  "Myeloid-state0","Myeloid-state5")
  
  # Correlation of states with genes
  state_cors <- cor(GetAssayData(visium_slide, assay = "SCT")[genes_to_cor, ] %>%
                      as.matrix() %>%
                      t(),
                    GetAssayData(visium_slide, assay = "cell_states") %>%
                      as.matrix() %>%
                      t()) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    dplyr::select_at(sel_states)
  
  colnames(state_cors) <- gsub(pattern = "-","_",colnames(state_cors))
  
  state_cors
  
}

# Take the median correlation between states and genes
# Later classify interactions as myof or regular and take the median again

gene_state_cors <- meta_data_slides[, c("spot_id", "orig.ident")] %>%
  group_by(orig.ident) %>%
  nest() %>%
  mutate() %>%
  mutate(file_path = paste0("./processed_visium/objects/", 
                            orig.ident,
                            ".rds")) %>%
  mutate(cor_res = map2(data, file_path, get_gene_state_cor)) %>%
  dplyr::select(orig.ident, cor_res) %>%
  unnest() %>%
  ungroup() %>%
  group_by(gene) %>%
  dplyr::select(-orig.ident) %>% 
  summarise_all(., median) %>% 
  pivot_longer(-gene) %>%
  dplyr::mutate(name = ifelse(name %in% c("Fib_state6","Myeloid_state5"),
                              "myof_int", "reg_int")) %>%
  dplyr::group_by(gene, name) %>%
  summarize(median_congruence = median(value)) %>%
  pivot_wider(names_from = name, values_from = median_congruence)

get_cor_df <-  function(dat, f_path) {
  
  visium_slide <- readRDS(f_path)[genes_to_cor, dat$spot_id]
  
  cor_markers <- cor(GetAssayData(visium_slide, assay = "SCT") %>%
                       as.matrix() %>%
                       t(), 
                     method = "spearman") %>%
    as.data.frame() %>%
    rownames_to_column("gene_a") %>%
    pivot_longer(-gene_a, names_to = "gene_b", values_to = "spearman_cor") %>%
    dplyr::filter(gene_a != gene_b) %>%
    left_join(all_candidates, by = c("gene_a" = "gene")) %>%
    dplyr::rename("state_a" = name) %>%
    left_join(all_candidates, by = c("gene_b" = "gene")) %>%
    dplyr::rename("state_b" = name)
  
}


all_cors <- meta_data_slides[, c("spot_id", "orig.ident")] %>%
  group_by(orig.ident) %>%
  nest() %>%
  mutate() %>%
  mutate(file_path = paste0("./processed_visium/objects/", 
                           orig.ident,
                           ".rds")) %>%
  mutate(cor_res = map2(data, file_path, get_cor_df))

candidate_interactions <- all_cors %>%
  dplyr::select(orig.ident, cor_res) %>%
  unnest(cor_res) %>%
  group_by(gene_a, gene_b, state_a, state_b) %>%
  summarize(mean_cor = mean(spearman_cor)) %>%
  ungroup() %>%
  group_by(state_a, state_b) %>%
  nest()


# We will correlate markers 

candidate_myof <- candidate_interactions %>%
  dplyr::filter(state_a == "m_state5",
                state_b == "f_state6") %>%
  unnest() %>%
  arrange(-mean_cor) %>%
  left_join(gene_state_cors[,c("gene", "myof_int")], 
            by = c("gene_a" = "gene")) %>%
  left_join(gene_state_cors[,c("gene", "myof_int")], 
            by = c("gene_b" = "gene")) %>%
  group_by(state_a, state_b, gene_a, gene_b) %>%
  dplyr::mutate(myof_score = median(c(myof_int.x,myof_int.y))) %>%
  dplyr::select(mean_cor, myof_score)

candidate_myof_filt <- candidate_myof %>%
  dplyr::filter((! gene_a %in% intersect(candidate_myof$gene_a, candidate_myof$gene_b)) &
                (! gene_b %in% intersect(candidate_myof$gene_a, candidate_myof$gene_b))) %>%
  dplyr::filter(mean_cor > 0.05) %>%
  arrange(-myof_score)

candidate_dormant <- candidate_interactions %>%
  dplyr::filter(state_a == "m_state0",
                state_b == "f_state5") %>%
  unnest() %>%
  arrange(-mean_cor) %>%
  left_join(gene_state_cors[,c("gene", "reg_int")], 
            by = c("gene_a" = "gene")) %>%
  left_join(gene_state_cors[,c("gene", "reg_int")], 
            by = c("gene_b" = "gene")) %>%
  group_by(state_a, state_b, gene_a, gene_b) %>%
  dplyr::mutate(reg_score = median(c(reg_int.x,reg_int.y))) %>%
  dplyr::select(mean_cor, reg_score)

candidate_dormant_filt <- candidate_dormant %>%
  dplyr::filter((! gene_a %in% intersect(candidate_dormant$gene_a, candidate_dormant$gene_b)) &
                  (! gene_b %in% intersect(candidate_dormant$gene_a, candidate_dormant$gene_b))) %>%
  dplyr::filter(mean_cor > 0.05) %>%
  arrange(-reg_score)

# Save tables

write.table(candidate_myof,
            file = "./results/rnascope/myeloid_fibro/IZ/candidate_myof.txt",
            col.names = T, row.names = F, quote = F, sep = "\t") 

write.table(candidate_myof_filt,
            file = "./results/rnascope/myeloid_fibro/IZ/candidate_myof_filt.txt",
            col.names = T, row.names = F, quote = F, sep = "\t") 
write.table(candidate_dormant,
            file = "./results/rnascope/myeloid_fibro/IZ/candidate_dormant.txt",
            col.names = T, row.names = F, quote = F, sep = "\t") 

write.table(candidate_dormant_filt, 
            file = "./results/rnascope/myeloid_fibro/IZ/candidate_dormant_filt.txt",
            col.names = T, row.names = F, quote = F, sep = "\t") 

## Sanity check
pdf("./results/rnascope/myeloid_fibro/IZ/sanity_check.pdf", height = 17, width = 10)

candidate_interactions %>%
  mutate(id = paste0(state_a,"_",state_b)) %>%
  unnest() %>%
  arrange(id, -mean_cor) %>%
  group_by(id) %>%
  dplyr::slice(1:50) %>%
  dplyr::filter(id %in% c("f_state5_f_state5",
                          "f_state6_f_state6",
                          "f_state5_f_state6",
                          "m_state5_m_state5",
                          "m_state0_m_state0",
                          "m_state5_m_state0",
                          "m_state5_f_state6",
                          "m_state5_f_state5",
                          "m_state0_f_state5",
                          "m_state0_f_state6")) %>%
  ggplot(.,aes(y = id, x = mean_cor)) +
  geom_boxplot()

dev.off()

# Check in some slides

param_df <- meta_data_slides[, c("spot_id", "orig.ident")] %>%
  group_by(orig.ident) %>%
  nest() %>%
  dplyr::select(-data) %>%
  mutate() %>%
  mutate(file_path = paste0("./processed_visium/objects/", 
                            orig.ident,
                            ".rds"))

library(cowplot)

plot_markers_myofib_list <- function(slide, f_path) {
  
  visium_slide <- readRDS(f_path)
  
  ids <- meta_data_slides[, c("spot_id", "orig.ident")] %>%
    dplyr::filter(orig.ident == slide) %>% 
    pull(spot_id)
  
  visium_slide <- visium_slide[,ids]
  
  myofib_list = map(seq(1:10), function(i) {
    
    ga <- candidate_myof_filt[i,"gene_a"][[1]]
    gb <- candidate_myof_filt[i,"gene_b"][[1]]
    
    SpatialFeaturePlot(visium_slide, features = c(ga, gb))
    
  })
  
  plot(cowplot::plot_grid(plotlist = myofib_list, ncol = 2))
  
}

plot_markers_dormant_list <- function(slide, f_path) {
  
  visium_slide <- readRDS(f_path)
  
  ids <- meta_data_slides[, c("spot_id", "orig.ident")] %>%
    dplyr::filter(orig.ident == slide) %>% 
    pull(spot_id)
  
  visium_slide <- visium_slide[,ids]
  
  myofib_list = map(seq(1:10), function(i) {
    
    ga <- candidate_dormant_filt[i,"gene_a"][[1]]
    gb <- candidate_dormant_filt[i,"gene_b"][[1]]
    
    SpatialFeaturePlot(visium_slide, features = c(ga, gb))
    
  })
  
  plot(cowplot::plot_grid(plotlist = myofib_list, ncol = 2))
  
}

pdf("./results/rnascope/myeloid_fibro/IZ/candidate_myof_filt.pdf", height = 17, width = 10)

walk2(param_df$orig.ident,param_df$file_path, plot_markers_myofib_list)

dev.off()

pdf("./results/rnascope/myeloid_fibro/IZ/candidate_dormant_filt.pdf", height = 17, width = 10)

walk2(param_df$orig.ident,param_df$file_path, plot_markers_dormant_list)

dev.off()
