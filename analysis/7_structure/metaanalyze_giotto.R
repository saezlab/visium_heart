# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' We look at the niche relationship that we get between groups

library(tidyverse)
library(survcomp)

# Loading patient grouping 

patient_meta <- read.table("./markers/visium_annotations_ext.txt",sep = '\t', header = T)
patient_order <- patient_meta$pid

# Loading neighborhood analysis 

giotto_results_folder <- "./results/niche_mapping/giotto/result_tables"

file_names <- list.files(giotto_results_folder, 
                         full.names = F)

sample_names <- gsub("_neighbor_res.rds", "", file_names)

giotto_res <- map(set_names(paste0(giotto_results_folder, "/",
                                  file_names),
                            sample_names), readRDS) %>%
  enframe("sample") %>%
  unnest()

# First we will take hetero interactions
# These interactions should be in most of tha samples
# These would describe the general structures in the heart

n_samples_prop = 28 * 0.6

general_interactions <- giotto_res %>%
  mutate(lowest_p = ifelse(p.adj_higher < p.adj_lower,
                           p.adj_higher, p.adj_lower)) %>%
  dplyr::filter(type_int == "hetero",
                lowest_p < 0.05) %>%
  group_by(unified_int) %>%
  mutate(n_appearances = length(lowest_p)) %>%
  dplyr::filter(n_appearances > n_samples_prop) %>%
  ungroup() %>%
  arrange(n_appearances, unified_int) %>%
  group_by(unified_int) %>%
  nest() %>%
  ungroup() %>%
  slice(1:30)

general_hmap <- general_interactions %>%
  unnest() %>%
  dplyr::mutate(pid = factor(pid,
                             levels = patient_order)) %>%
  ggplot(aes(x = pid, 
             y = unified_int, 
             fill = PI_value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(midpoint = 0, na.value = 'white')

pdf("./results/niche_mapping/giotto/plots/compiled_interactions_general.pdf", height = 5)

plot(general_hmap)

dev.off()

# Now, each patient group defined by pseudobulk may have a particular set of interactions that are useful
# Assumption: PI_value = 0 when the test wasn't done

run_wilcox_all <- function(prop_data) {
  
  prop_data_group <- prop_data[["patient_group"]] %>%
    unique() %>%
    set_names()
  
  map(prop_data_group, function(g) {
    
    test_data <- prop_data %>%
      mutate(test_group = ifelse(patient_group == g,
                                 "target", "rest"))
    
    wilcox.test(PI_value ~ test_group, 
                data = test_data,
                alternative = "two.sided") %>%
      broom::tidy()
  }) %>% enframe("patient_group") %>%
    unnest()
  
}

wilcoxon_test_niche_interactions <- giotto_res %>%
  left_join(patient_meta, by = c("sample" = "sample_id")) %>%
  dplyr::select(pid, unified_int, PI_value) %>%
  pivot_wider(names_from = pid, values_from = PI_value, values_fill = 0) %>%
  pivot_longer(-unified_int, names_to = "pid", values_to = "PI_value") %>%
  left_join(patient_meta) %>%
  group_by(unified_int) %>%
  nest() %>%
  mutate(wilcox_res = map(data, run_wilcox_all)) %>%
  dplyr::select(unified_int, wilcox_res) %>%
  unnest() %>%
  mutate(corr_pval = p.adjust(p.value)) %>%
  dplyr::filter(corr_pval <= 0.15) %>%
  arrange(patient_group, corr_pval)

# Get hetero interactions of these tests

group_interactions <- giotto_res %>%
  dplyr::filter(type_int == "hetero") %>%
  left_join(patient_meta, by = c("sample" = "sample_id")) %>%
  dplyr::select(pid, unified_int, PI_value) %>%
  pivot_wider(names_from = pid, values_from = PI_value, values_fill = 0) %>%
  pivot_longer(-unified_int, names_to = "pid", values_to = "PI_value") %>%
  left_join(patient_meta) %>%
  group_by(unified_int) %>%
  dplyr::filter(unified_int %in% wilcoxon_test_niche_interactions$unified_int) %>%
  nest()


pdf("./results/niche_mapping/giotto/plots/interaction_contrasts.pdf", height = 4, width = 4)

walk2(group_interactions$unified_int, group_interactions$data, function(interaction, int_data) {
  
  plot(int_data %>%
    ggplot(aes(x = patient_group, y = PI_value , fill = patient_group)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90)) +
    ylab("PI_value") +
    ggtitle(interaction))
  
})

dev.off()

  mutate(lowest_p = ifelse(p.adj_higher < p.adj_lower,
                           p.adj_higher, p.adj_lower)) %>%
  dplyr::filter(type_int == "hetero",
                lowest_p < 0.05) %>%
  group_by(unified_int) %>%
  mutate(n_appearances = length(lowest_p)) %>%
  dplyr::filter(n_appearances > n_samples_prop) %>%
  ungroup() %>%
  arrange(n_appearances, unified_int) %>%
  group_by(unified_int) %>%
  nest() %>%
  ungroup() %>%
  slice(1:30)



# Visualizing hetero relationships (higher)

pos_dependencies <- giotto_res %>%
  dplyr::filter(type_int == "hetero") %>%
  dplyr::select(unified_int, sample, PI_value, p.adj_higher) %>%
  left_join(patient_meta, by = c("sample" = "sample_id")) %>%
  dplyr::mutate(pid = factor(pid,
                                levels = patient_order))

interaction_order <- pos_dependencies %>%
  group_by(unified_int) %>%
  summarise(n_occur = length(PI_value),
            mean_int_score = mean(PI_value)) %>%
  dplyr::filter(n_occur >= 5) %>%
  arrange(-mean_int_score) %>%
  slice(1:30) %>%
  pull(unified_int)

int_hmap <- pos_dependencies %>%
  dplyr::filter(unified_int %in% interaction_order) %>%
  dplyr::select(pid, unified_int, PI_value) %>%
  pivot_wider(names_from = pid, values_from = PI_value) %>%
  pivot_longer(-unified_int, names_to = "pid", values_to = "PI_value") %>%
  dplyr::mutate(pid = factor(pid,
                             levels = patient_order),
                unified_int = factor(unified_int,
                                     levels = interaction_order)) %>%
  ggplot(aes(x = pid, 
             y = unified_int, 
             fill = PI_value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(midpoint = 0, na.value = 'white')

pdf("./results/niche_mapping/giotto/plots/compiled_interactions_hmap.pdf", height = 8)

plot(int_hmap)

dev.off()

# Visualizing lower

neg_dependencies <- giotto_res %>%
  dplyr::filter(type_int == "hetero") %>%
  dplyr::select(unified_int, sample, PI_value, p.adj_lower) %>%
  left_join(patient_meta, by = c("sample" = "sample_id")) %>%
  dplyr::mutate(pid = factor(pid,
                             levels = patient_order))

interaction_order <- neg_dependencies %>%
  group_by(unified_int) %>%
  summarise(n_occur = length(PI_value),
            mean_int_score = mean(PI_value)) %>%
  dplyr::filter(n_occur >= 5) %>%
  arrange(mean_int_score) %>%
  slice(1:30) %>%
  pull(unified_int)

int_hmap <- neg_dependencies %>%
  dplyr::filter(unified_int %in% interaction_order) %>%
  dplyr::select(pid, unified_int, PI_value) %>%
  pivot_wider(names_from = pid, values_from = PI_value) %>%
  pivot_longer(-unified_int, names_to = "pid", values_to = "PI_value") %>%
  dplyr::mutate(pid = factor(pid,
                             levels = patient_order),
                unified_int = factor(unified_int,
                                     levels = interaction_order)) %>%
  ggplot(aes(x = pid, 
             y = unified_int, 
             fill = PI_value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(midpoint = 0, na.value = 'white')

pdf("./results/niche_mapping/giotto/plots/compiled_interactions_hmap_lower.pdf", height = 6)

plot(int_hmap)

dev.off()

# Try alternative wth complexhmap

int_mat <- pos_dependencies %>%
  dplyr::filter(unified_int %in% interaction_order) %>%
  dplyr::select(pid, unified_int, PI_value) %>%
  pivot_wider(names_from = pid, values_from = PI_value) %>%
  mutate(unified_int = as.character(unified_int) %>%
           gsub("--","_",.)) %>%
  column_to_rownames("unified_int") %>%
  as.matrix()

int_mat[is.na(int_mat)] <- 0

pdf("./results/niche_mapping/giotto/plots/compiled_interactions_hmap_alt.pdf", height = 6)
(ComplexHeatmap::Heatmap(int_mat[, patient_order],
                             cluster_columns = FALSE))
dev.off()

# You can try to do this per group

group_meta_res <- pos_dependencies %>% 
  group_by(patient_group) %>%
  nest() %>%
  mutate(filt_data = map(data, function(dat) {
    
    filt_dat <- dat %>% group_by(unified_int) %>%
      mutate(n_occur = length(PI_value),
                mean_int_score = mean(PI_value)) %>%
      dplyr::filter(n_occur >= 3) %>%
      summarise(meta_pval = median(p.adj_higher))
  
    
  })) %>% 
  dplyr::select(patient_group, filt_data) %>%
  unnest()


pdf("./results/niche_mapping/giotto/plots/compiled_interactions_hmap_median_pval.pdf", height = 4, width = 4)
plot(group_meta_res %>% 
  dplyr:::filter(meta_pval < 0.15) %>%
  dplyr::mutate(meta_pval = ifelse(meta_pval < 0.00001, 0.00001, meta_pval)) %>%
  ggplot(aes(x = patient_group, y = unified_int, fill = -log10(meta_pval))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
dev.off()


