# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
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

# Now let's get all the networks

regulons <- tibble(grn_name = list.files("./reg_nets/processed/")) %>%
  mutate(cell_type = gsub("[.]txt", "", grn_name)) %>%
  mutate(grn = map(paste0("./reg_nets/processed/", grn_name), read_table2)) %>%
  dplyr::select(-grn_name) %>%
  unnest() %>%
  dplyr::mutate(source = paste0(cell_type, 
                                "_",
                                source)) %>%
  dplyr::select(source, target, likelihood, mor)

min_targets <- 5

regulons <- regulons %>% 
  dplyr::filter(target %in% rownames(norm_pseudobulk_data_counts))

filtered_regulons <- regulons %>% 
  dplyr::select(source, target) %>%
  group_by(source) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  mutate(n_targets = map_dbl(data, length)) %>%
  dplyr::filter(n_targets >= min_targets) %>%
  dplyr::select(-n_targets) %>%
  deframe() %>%
  names()

regulons <- regulons %>%
  dplyr::filter(source %in% filtered_regulons)

# Generate matrix of t-values

tf_acts <- decoupleR::run_wmean(mat = norm_pseudobulk_data_counts, 
                                     network = regulons,
                                     .source = "source",
                                     .target = "target",
                                     .mor = "mor",
                                     .likelihood = "likelihood",
                                     times = 1000)


# Use weighted mean (normalized)

tf_acts <- tf_acts %>%
  dplyr::filter(statistic == "norm_wmean") %>%
  dplyr::mutate(ct_TF = strsplit(source, "_") %>%
                  map_chr(., ~ .x[[1]])) %>%
  dplyr::mutate(source = strsplit(source, "_") %>%
                  map_chr(., ~ .x[[2]]))

# Filter TFs by correspondant TF
# Nest by ct (condition)
# if condition in ct_TF then use that value
# else, use the mean

tf_acts <- tf_acts %>%
  dplyr::select(-statistic) %>%
  dplyr::group_by(condition) %>%
  nest() %>%
  mutate(tf_acts = map2(condition, data, function(ct, dat) {
    
    cts_dat <- dat[["ct_TF"]]
    
    if(ct %in% cts_dat) {
      
      tf_info <- dat %>%
        dplyr::filter(ct_TF == ct) %>%
        dplyr::select(-c("ct_TF"))
        
      return(tf_info)
    } else {
      
      tf_info <- dat %>%
        group_by(source) %>%
        summarise(score = mean(score),
                  p_value = mean(p_value)) 
      
    }

  }))

# Make a matrix to be able to correlate
tf_acts <- tf_acts %>%
  dplyr::select(tf_acts) %>%
  unnest() %>%
  ungroup()

tf_acts_mat <- tf_acts %>%
  dplyr::select(-p_value) %>%
  pivot_wider(values_from = score,
              names_from = condition,
              values_fill = 0) %>%
  column_to_rownames("source") %>%
  as.matrix()


# Load HINT results --------------------

hint_res <- read_csv("./HINT/atlas/tf_activity.csv") %>%
  dplyr::rename("source" = X1) %>%
  dplyr::rename("PC" = Pericyte) %>%
  dplyr::mutate(source = strsplit(source, "[.]") %>%
                  map_chr(., ~.x[[3]]))

hint_res <- hint_res %>%
  pivot_longer(-source, 
               names_to = "condition",
               values_to = "score_hint") %>%
  left_join(tf_acts %>%
              dplyr::select(-p_value),
            by = c("condition", "source")) %>%
  dplyr::mutate(score = ifelse(is.na(score), 0, score)) %>%
  group_by(source) %>%
  mutate(std_score = (score - mean(score))/sd(score)) %>%
  dplyr::filter(!is.nan(std_score)) %>%
  dplyr::mutate(source = toupper(source))

# Here we make all genes upper
TF_annotation <- hint_res %>%
  dplyr::filter(score_hint == max(score_hint)) %>%
  dplyr::select(source, condition)

# Save HINT/Regulome activities

write_csv(hint_res, "./results/atacvrna/HINT_regulome_unif_acts.csv")

# Let's do some summaries

ct_correlation <- hint_res %>%
  group_by(condition) %>%
  nest() %>%
  mutate(TFcor = map(data, function(dat) {
    
    cor.test(dat$score_hint, dat$std_score, method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  select(TFcor) %>%
  unnest()

write_csv(ct_correlation, "./results/atacvrna/HINT_regulome_ct_cors.csv")

tf_correlation <- hint_res %>%
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
  
write_csv(tf_correlation, "./results/atacvrna/HINT_regulome_tf_cors.csv")

# Make a selection and create heatmaps
selected_TFs <- tf_correlation %>%
  group_by(condition) %>%
  dplyr::slice(1:5) %>%
  pull(source)

plt_df <- hint_res %>%
  dplyr::filter(source %in% selected_TFs) %>%
  dplyr::mutate(source = factor(source,
                                levels = selected_TFs))

write_csv(plt_df, file = "./results/atacvrna/atlasTFconsistency.csv")

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

fig_panel <- cowplot::plot_grid(hmapA, hmapB, nrow = 1)

pdf("./results/atacvrna/atlasTFconsistency.pdf", height = 6, width = 8)
plot(fig_panel)
dev.off()

# Now let's check cell correlations

ct_correlation_filt <- hint_res %>%
  dplyr::filter(source %in% unique(tf_correlation$source)) %>%
  group_by(condition) %>%
  nest() %>%
  mutate(TFcor = map(data, function(dat) {
    
    cor.test(dat$score_hint, dat$std_score, method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  select(TFcor) %>%
  unnest()

write_csv(ct_correlation_filt, "./results/atacvrna/HINT_regulome_ct_cors_filt.csv")

# Let's try a second test, where we include gene expression

scaled_expression <- t(norm_pseudobulk_data_counts[, hint_res$condition %>% 
                                                     unique()]) %>% 
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("condition") %>%
  pivot_longer(-condition, names_to = "source", values_to = "expr")

all_mods <- hint_res %>%
  left_join(scaled_expression, by = c("source", "condition")) %>%
  dplyr::filter(!is.na(expr))

ct_correlation_all <- all_mods %>%
  #dplyr::filter(source %in% unique(tf_correlation$source)) %>%
  group_by(condition) %>%
  nest() %>%
  mutate(hint_regulon_cor = map(data, function(dat) {
    
    cor.test(abs(dat$score_hint), abs(dat$std_score), method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  mutate(hint_expr_cor = map(data, function(dat) {
    
    cor.test(abs(dat$score_hint), abs(dat$expr), method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  mutate(regulon_expr_cor = map(data, function(dat) {
    
    cor.test(abs(dat$std_score), abs(dat$expr), method = "spearman") %>%
      broom::tidy()
    
  }))

ct_correlation_all_sign <- all_mods %>%
  dplyr::filter(source %in% unique(tf_correlation$source)) %>%
  group_by(condition) %>%
  nest() %>%
  mutate(hint_regulon_cor = map(data, function(dat) {
    
    cor.test((dat$score_hint), (dat$std_score), method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  mutate(hint_expr_cor = map(data, function(dat) {
    
    cor.test((dat$score_hint), (dat$expr), method = "spearman") %>%
      broom::tidy()
    
  })) %>%
  mutate(regulon_expr_cor = map(data, function(dat) {
    
    cor.test((dat$std_score), (dat$expr), method = "spearman") %>%
      broom::tidy()
    
  }))

saveRDS(ct_correlation_all_sign, "./results/atacvrna/TF_modality_consistency.rds")

ct_plts <- ct_correlation_all %>%
  dplyr::select(-data) %>%
  pivot_longer(-condition) %>%
  unnest() %>%
  ggplot(aes(x = name, y = condition, fill = estimate)) +
  geom_tile() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)) +
  scale_fill_gradient2() +
  ylab("Cell type") + xlab("Modality comparison")

ct_plts_sign <- ct_correlation_all_sign %>%
  dplyr::select(-data) %>%
  pivot_longer(-condition) %>%
  unnest() %>%
  ggplot(aes(x = name, y = condition, fill = estimate)) +
  geom_tile() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)) +
  scale_fill_gradient2() +
  ylab("Cell type") + xlab("Modality comparison")


pdf("./results/atacvrna/TF_modality_consistency.pdf", height = 4, width = 3)

plot(ct_plts)
plot(ct_plts_sign)

dev.off()

