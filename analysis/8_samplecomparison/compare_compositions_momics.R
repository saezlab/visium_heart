# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we compare the compositions of each patient group
#' using each modality or all of them

library(tidyverse)

# Read individual compositions

atac_props <- read.csv("./results/compositions/atac_compositions.txt", sep = "\t")
rna_props <- read.csv("./results/compositions/snrna_compositions.txt", sep = "\t")
spatial_props <- read.csv("./results/compositions/spatial_compositions.txt", sep = "\t")

# patient grouping
pat_anns <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  dplyr::select(patient_id, patient_group) %>%
  unique()

# All proportions
mi_props <- left_join(spatial_props, rna_props, 
                      by = c("patient_id", "cell_type")) %>%
  left_join(atac_props,
            by = c("patient_id", "cell_type")) %>%
  mutate(sn_n_cells = ifelse(is.na(sn_n_cells), 0, sn_n_cells),
         sn_prop_cells = ifelse(is.na(sn_prop_cells), 0, sn_prop_cells),
         atac_n_cells = ifelse(is.na(atac_n_cells), 0, atac_n_cells),
         atac_prop_cells = ifelse(is.na(atac_prop_cells), 0, atac_prop_cells)) %>%
  left_join(pat_anns)

# Correlate the mean proportions per group

mi_props_mean <- mi_props %>%
  group_by(patient_group, cell_type) %>%
  summarize(mean_prop_rna = mean(sn_prop_cells),
            mean_prop_spatial = mean(sp_prop_cells),
            mean_prop_atac = mean(atac_prop_cells)) %>%
  nest() %>%
  mutate(cor_res = map(data, function(dat) {
    
    correlation = bind_rows(cor.test(log1p(dat$mean_prop_rna), 
                                     log1p(dat$mean_prop_atac),
                                     method = "spearman") %>%
                              broom::tidy(),
                            cor.test(log1p(dat$mean_prop_rna), 
                                     log1p(dat$mean_prop_spatial),
                                     method = "spearman") %>%
                              broom::tidy(),
                            cor.test(log1p(dat$mean_prop_atac), 
                                     log1p(dat$mean_prop_spatial),
                                     method = "spearman") %>%
                              broom::tidy())
    
    bind_cols(comparison = c("rna_atac_cor", "rna_spatial_cor", "atac_spatial_cor"),
              correlation = correlation)
    
  })) %>%
  select(cor_res) %>%
  unnest() 

write.table(mi_props_mean, file = "./results/compositions/patientgroup_compositions_cor.txt", 
            col.names = T, row.names = F, quote = F, sep = "\t")



pdf(file = "./results/compositions/patientgroup_compositions_cor.pdf", height = 3, width = 3.5)

mi_props_mean %>%
  ggplot(aes(x = patient_group, y = comparison, fill = estimate)) +
  geom_tile() +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  scale_fill_gradient(low = "#ffd89b", high = "#19547b", limits = c(0,1)) 


dev.off()


# Let's create 4 panels

pdf(file = "./results/sample_comparison/compositions/rna_cellcomps_comparisons.pdf", height = 4, width = 5)

rna_compositions_plt <- ggplot(mi_props, aes(x = patient_group, y = sn_prop_cells, color = patient_group)) +
  geom_boxplot() +
  facet_wrap(.~cell_type, ncol = 4, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("snRNA-seq estimated compositions") +
  xlab("")

plot(rna_compositions_plt)

dev.off()

pdf(file = "./results/sample_comparison/compositions/atac_cellcomps_comparisons.pdf", height = 4, width = 5)

atac_compositions_plt <- ggplot(mi_props, aes(x = patient_group, y = atac_prop_cells, color = patient_group)) +
  geom_boxplot() +
  facet_wrap(.~cell_type, ncol = 4, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("snATAC-seq estimated compositions") +
  xlab("")

plot(atac_compositions_plt)

dev.off()

pdf(file = "./results/sample_comparison/compositions/spatial_cellcomps_comparisons.pdf", height = 4, width = 5)

spatial_compositions_plt <- ggplot(mi_props, aes(x = patient_group, y = sp_prop_cells, color = patient_group)) +
  geom_boxplot() +
  facet_wrap(.~cell_type, ncol = 4, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("deconvolution estimated compositions") +
  xlab("")

plot(spatial_compositions_plt)

dev.off()

# Now take the mean per patient

multiview_props <- mi_props %>%
  pivot_longer(-c("patient_id", "cell_type", "patient_group")) %>%
  dplyr::filter(grepl("prop", name)) %>%
  group_by(patient_id, cell_type) %>%
  summarise(multiview_mean_prop = mean(value)) %>%
  left_join(pat_anns)

multiview_props %>% write_csv("./results/sample_comparison/compositions/momics_cellcomps.csv")

# Statistical tests:

kwallis_comps <- multiview_props %>% 
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = map(data, function(dat){
    
    kruskal.test(multiview_mean_prop ~ patient_group, data = dat) %>%
      broom::tidy()
    
  })) %>%
  select(wres) %>%
  unnest() %>%
  ungroup() %>%
  mutate(corr_p = p.adjust(p.value))


write.table(kwallis_comps,
              file = "./results/sample_comparison/compositions/kruskall_wallis_momics_cellcomps.txt", 
              col.names = T, row.names = F, quote = F, sep = "\t")


winners <- kwallis_comps %>%
  dplyr::filter(corr_p <= 0.1) %>%
  pull(cell_type)

pdf(file = "./results/sample_comparison/compositions/momics_cellcomps.pdf", height = 3, width = 5)

multiview_compositions_plt <- multiview_props %>%
  dplyr::filter(cell_type %in% winners) %>%
  ggplot(., aes(x = patient_group, y = multiview_mean_prop, color = patient_group)) +
  geom_boxplot() +
  facet_wrap(.~cell_type, ncol = 4, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  ylab("multiomics estimated\n compositions") +
  xlab("")

plot(multiview_compositions_plt)

dev.off()

# Generate the clustering based on cell compositions

multiview_props_mat <- multiview_props %>%
  dplyr::select(-patient_group) %>%
  pivot_wider(names_from = cell_type, values_from = multiview_mean_prop) %>%
  as.data.frame() %>%
  column_to_rownames("patient_id") %>%
  as.matrix()

# Generate ILR transformation
baseILR <- ilrBase(x = multiview_props_mat,
                   method = "basic")

cell_ilr <- as.matrix(ilr(multiview_props_mat, baseILR))

colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

gex_hclust <- eclust(cell_ilr, "hclust", k = 3)

# Make color palette

color_palette <- tibble(patient_id = gex_hclust$labels[gex_hclust$order]) %>%
  left_join(pat_anns[,c("patient_group", "patient_id")] %>% unique()) %>%
  left_join(tibble(patient_group = c("group_1", "group_2", "group_3"),
                   col = c("red", "darkgreen", "blue")))

pdf("./results/sample_comparison/compositions/momics_cellcomps_patient_clustering.pdf", height = 6, width = 5)

plot(fviz_dend(gex_hclust, 
               rect = TRUE, 
               label_cols = color_palette$col,
               k_colors = rep("black",3)))

dev.off()


