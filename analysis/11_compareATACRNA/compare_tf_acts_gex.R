# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate TF activities on the pseudobulk data 
#' using corrected wsum and 

library(scater)
library(tidyverse)
library(decoupleR)
library(cowplot)
source("./analysis/utils/pseudobulk_utils.R")

# Reading the pseudobulk data of the atlas

pseudobulk_data <- readRDS("./processed_snrnaseq/integration/ps_integrated_rnasamples_ann.rds")[[1]][["gex"]]
pseudobulk_data_counts <- assay(pseudobulk_data)
colnames(pseudobulk_data_counts) <- pseudobulk_data$slide.meta.data...vars.

# Normalization and filtering parameters... Are these too strict?

pseudobulk_data_counts <- edgeR_filtering(pseudobulk_data_counts, 
                                          min.count = 10,
                                          min.prop = 0.05,
                                          min.total.count = 10)


norm_pseudobulk_data_counts <- cpm_norm(pseudobulk_data_counts)

#scaled_pseudobulk_data_counts <- scale(t(norm_pseudobulk_data_counts)) %>%
#  t()

# Reading TF acts

hint_res <- read_csv("./results/atacvrna/HINT_regulome_unif_acts.csv")

TF_annotation <- hint_res %>%
  group_by(source) %>%
  dplyr::filter(score_hint == max(score_hint)) %>%
  dplyr::select(source, condition)

# Filtering by gene expression

scaled_expression <- t(norm_pseudobulk_data_counts[, hint_res$condition %>% 
                                                     unique()]) %>% 
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("condition") %>%
  pivot_longer(-condition, names_to = "source", values_to = "expr")

all_mods <- hint_res %>%
  left_join(scaled_expression, by = c("source", "condition")) %>%
  dplyr::filter(!is.na(expr))


tf_correlation <- all_mods %>%
  group_by(source) %>%
  nest() %>%
  mutate(TFcor = map(data, function(dat) {
    
    cor.test(dat$score_hint, dat$std_score, method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  select(TFcor) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::filter(estimate > 0.3) %>%
  left_join(TF_annotation) %>%
  dplyr::arrange(condition, -estimate)

# Make a selection and create heatmaps
selected_TFs <- tf_correlation %>%
  group_by(condition) %>%
  dplyr::slice(1:5) %>%
  pull(source)

plt_df <- all_mods %>%
  dplyr::filter(source %in% selected_TFs) %>%
  dplyr::mutate(source = factor(source,
                                levels = selected_TFs))

write_csv(plt_df, file = "./results/atacvrna/atlasTFconsistency_v2.csv")

hmapA <- plt_df %>%
  ggplot(aes(x = condition, y = source, fill = score_hint)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  scale_fill_gradient2(limits = c(-5,5)) +
  ylab("TF (HINT activities)") + xlab("")

hmapB <- plt_df %>%
  ggplot(aes(x = condition, y = source, fill = std_score)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  scale_fill_gradient2(limits = c(-5,5)) +
  ylab("TF (Regulon mean expression)") + xlab("")

hmapC <- plt_df %>%
  ggplot(aes(x = condition, y = source, fill = expr)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  scale_fill_gradient2(limits = c(-5,5)) +
  ylab("TF expression") + xlab("")

fig_panel <- cowplot::plot_grid(hmapA, hmapB, hmapC, nrow = 1)

pdf("./results/atacvrna/atlasTFconsistency_v2.pdf", height = 6, width = 12)
plot(fig_panel)
dev.off()


