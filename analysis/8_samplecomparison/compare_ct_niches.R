# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' This calculates the differences of niche composition between patient groups
library(tidyverse)
library(Seurat)

cluster_counts <- read_csv(file = "./results/niche_mapping/ct_niches/niche_prop_summary.csv")
niche_summary_pat <- read_table2(file = "./results/niche_mapping/ct_niches/niche_summary_pat.txt")
cluster_info <- readRDS("./results/niche_mapping/ct_niches/niche_annotation_ct.rds")
patient_info <- readRDS("./markers/visium_patient_anns_revisions.rds")

cluster_info <- cluster_info %>%
  dplyr::mutate(sample_id = strsplit(row_id, "[..]") %>%
                  map_chr(., ~ .x[[1]])) %>%
  dplyr::select(-niche) %>%
  left_join(patient_info)
  

patient_group_inf <- patient_info %>%
  dplyr::select(patient_id, patient_group) %>%
  unique()

# First get the descriptors of the niches ----------------------------------------------------

niche_summary <- niche_summary_pat %>%
  ungroup() %>%
  group_by(ct_niche, cell_type) %>%
  summarise(patient_median_ct_prop = median(median_ct_prop))

# Data manipulation to have clustered data

niche_summary_mat <- niche_summary %>%
  pivot_wider(values_from = patient_median_ct_prop, 
              names_from =  cell_type, values_fill = 0) %>%
  column_to_rownames("ct_niche") %>%
  as.matrix()

niche_order <- hclust(dist(niche_summary_mat))
niche_order <- niche_order$labels[niche_order$order]

ct_order <- hclust(dist(t(niche_summary_mat)))
ct_order <- ct_order$labels[ct_order$order]

# Find characteristic cell types of each niche
# We have per patient the proportion of each cell-type in each niche

run_wilcox_up <- function(prop_data) {
  
  prop_data_group <- prop_data[["ct_niche"]] %>%
    unique() %>%
    set_names()
  
  map(prop_data_group, function(g) {
    
    test_data <- prop_data %>%
      mutate(test_group = ifelse(ct_niche == g,
                                 "target", "rest")) %>%
      mutate(test_group = factor(test_group,
                                 levels = c("target", "rest")))
    
    wilcox.test(median_ct_prop ~ test_group, 
                data = test_data,
                alternative = "greater") %>%
      broom::tidy()
  }) %>% enframe("ct_niche") %>%
    unnest()
  
}

wilcoxon_res <- niche_summary_pat %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = map(data, run_wilcox_up)) %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value)) %>%
  mutate(significant = ifelse(p_corr <= 0.15, "*", ""))

# Plot 

mean_ct_prop_plt <- niche_summary %>%
  left_join(wilcoxon_res, by = c("ct_niche", "cell_type")) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order),
         ct_niche = factor(ct_niche, levels = niche_order)) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
  ungroup() %>%
  ggplot(aes(x = cell_type, y = ct_niche, fill = scaled_pat_median)) +
  geom_tile() +
  geom_text(aes(label = significant)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.y = element_text(size=12)) +
  scale_fill_gradient(high = "#ffd89b", low = "#19547b") +
  ggplot2::coord_equal() +
  ylab("")

pdf("./results/niche_mapping/ct_niches/ct_niches_description.pdf", height = 4, width = 5)

plot(mean_ct_prop_plt)

niche_summary %>%
  left_join(wilcoxon_res, by = c("ct_niche", "cell_type")) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order),
         ct_niche = factor(ct_niche, levels = niche_order)) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
  ungroup() %>%
  write_csv(., "./results/niche_mapping/ct_niches/ct_niches_description.csv")

dev.off()

# Calculate compositions of niches per sample -------------------------------------

niche_info <- cluster_info %>%
  dplyr::group_by(patient_id, ct_niche) %>%
  summarize(ncells = length(ct_niche)) %>% # spots per patient of each group
  dplyr::mutate(all_sample_cells = sum(ncells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cell_prop = ncells/all_sample_cells)

# Filter slides whenever a niche represents more than 80% of spots

high_qc_pats <- niche_info %>%
  group_by(patient_id) %>%
  summarize(max_comp =  max(cell_prop)) %>%
  dplyr::filter(max_comp < 0.8) %>%
  pull(patient_id)

complete_niche_info <- niche_info %>%
  dplyr::filter(patient_id %in% high_qc_pats) %>%
  dplyr::select(patient_id, ct_niche, cell_prop) %>%
  tidyr::complete(patient_id, ct_niche, fill = list("cell_prop" = 0)) %>%
  left_join(dplyr::select(patient_info, patient_id, patient_group) %>% unique())

# Make group comparisons of compositions ----------------------------------------------

niche_test_kw <- complete_niche_info %>%
  dplyr::select(ct_niche, cell_prop, patient_group) %>%
  group_by(ct_niche) %>%
  nest() %>%
  mutate(kw_res = map(data, function(dat) {
    
    kruskal.test(cell_prop ~ patient_group, 
                 data = dat) %>%
      broom::tidy()
    
    
  })) %>%
  dplyr::select(ct_niche, kw_res) %>%
  unnest() %>%
  ungroup() %>%
  mutate(corr_pval = p.adjust(p.value))

write.table(niche_test_kw, 
            file = "./results/sample_comparison/niche/kw_ct_niches_patients.txt", 
            col.names = T, 
            row.names = F, 
            quote = F, 
            sep = "\t")

selected_niches <- niche_test_kw %>%
  dplyr::filter(corr_pval <= .1) %>%
  pull(ct_niche)

pdf("./results/sample_comparison/niche/kw_ct_niches_patients.pdf", height = 3, width = 5)

complete_niche_info %>%
  dplyr::filter(ct_niche %in% selected_niches) %>%
  ggplot(aes(x = patient_group, y = cell_prop, fill = patient_group)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 90)) +
  facet_wrap(.~ ct_niche, nrow = 2, scales = "free_y") +
  ylab("niche proportion")

dev.off()

# Rest of the figures

# First, give an example of a niche

selected_slides <- complete_niche_info %>%
  dplyr::filter(ct_niche %in% c("niche_4", "niche_9")) %>%
  arrange(ct_niche, -cell_prop) %>%
  group_by(ct_niche) %>%
  slice(1:6) %>%
  left_join(patient_info)

# Selected slice:

pdf("./results/sample_comparison/niche/niche_4_example.pdf")

walk(selected_slides$sample_id, function(slide_f) {
  
  print(slide_f)
  
  visium_slide <- readRDS(paste0("./processed_visium/objects/", slide_f, ".rds"))
  
  DefaultAssay(visium_slide) <- "c2l_props"
  plot(SpatialFeaturePlot(visium_slide, features = c("Myeloid", "Mast", "Fib")))
  
  visium_slide <- visium_slide[, visium_slide$opt_clust_integrated == "niche_4"]
  plot(SpatialDimPlot(visium_slide,group.by = "opt_clust_integrated")  +
         scale_fill_manual(values = "#97BC62"))
  
})

dev.off()


pdf("./results/sample_comparison/niche/niche_density_example.pdf")

walk(selected_slides$sample_id, function(slide_f) {
  
  print(slide_f)
  
  visium_slide <- readRDS(paste0("./processed_visium/objects/", slide_f, ".rds"))
  
  DefaultAssay(visium_slide) <- "c2l_props"
  
  visium_slide <- visium_slide[, visium_slide$opt_clust_integrated %in% c("niche_9", "niche_4")]
  
  plot(SpatialDimPlot(visium_slide,group.by = "opt_clust_integrated") +
         scale_fill_manual(values = c("#97BC62FF", "#2C5F2D")))
  
})

dev.off()



