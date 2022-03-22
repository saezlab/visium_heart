# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will plot the death scores with UMI numbers

library(Seurat)
library(tidyverse)
library(viridis)
library(ggpubr)
source("./analysis/utils/funcomics.R")

snrna_qc <- read_csv("./processed_snrnaseq/initial_qc/all_qcs_after_filtering_snrna.csv")

pat_anns <- read_csv("markers/visium_patient_anns_revisions.csv") %>%
  left_join(snrna_qc, by = c("patient_region_id" = "patient_id"))

cell_death_scores <- readRDS(file = "./results/cell_death/cell_death_scores.rds")

sets <- c("BIOCARTA-DEATH-PATHWAY", 
          "REACTOME-REGULATED-NECROSIS") %>%
  set_names()


cell_death_scores  <- cell_death_scores %>% 
  mutate(death_scores = map(death_scores, function(dat) {
    dat[sets, ]
    }))

spatial_death <- function(sample_id, slide_file, death_scores) {
  
  print(sample_id)
  
  visium_slide <- readRDS(slide_file)
  
  for(gset in sets) {
    visium_slide[[gset]] <- death_scores[gset, colnames(visium_slide)]
  }
  
  scores <- visium_slide@meta.data[,c("nCount_Spatial", make.names(sets))] %>%
    rownames_to_column("spot_id")
  
  cor_df <- scores %>%
    pivot_longer(-c("nCount_Spatial", "spot_id")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(cor_res = map(data, function(dat) {
      cor.test(dat$nCount_Spatial, dat$value, method = "spearman") %>%
        broom::tidy()
    })) %>%
    dplyr::select(cor_res) %>%
    unnest() %>%
    mutate(sid = slide_file)
  
  pdf(paste0("./results/cell_death/sptl_plts/", sample_id, ".pdf"), height = 6, width = 6)
  
  all_plots <- SpatialPlot(visium_slide, 
                           features = c("nCount_Spatial", 
                                        make.names(sets)),
                           ncol = 2,
                           max.cutoff = "q99",combine = F)
  
  all_plots = map(all_plots, ~ .x + scale_fill_viridis(option = "A"))
  
  plot(cowplot::plot_grid(plotlist = all_plots, ncol = 2, align = "hv"))
  

  dev.off()
  
  return(cor_df)
  
}

# Correlation between UMIs and death scores

cor_res <- pmap(cell_death_scores, spatial_death)
names(cor_res) = cell_death_scores$sample_id

cor_res %>% enframe(name = "sample_id") %>% unnest() %>% 
  #dplyr::filter(name == "REACTOME.REGULATED.NECROSIS") %>%
  arrange(name, estimate) %>%
  left_join(pat_anns) %>%
  write_csv("./results/cell_death/cor_UMI_death.csv")


# Calculate mean scores

mean_scores <- map(set_names(cell_death_scores$death_scores, 
                             cell_death_scores$sample_id), function(x) {
                               x[is.infinite(x)] <- NA
                               m_s <- rowMeans(x, na.rm = T)
                               tibble(gset = names(m_s), score = m_s)
                             }) %>%
  enframe() %>%
  unnest() %>%
  left_join(pat_anns %>% dplyr::select(sample_id, major_labl) %>% unique(),
            by = c("name" = "sample_id")) %>%
  group_by(gset) %>%
  nest()

unnest(mean_scores) %>%
  write_csv("./results/cell_death/wmean_score_all.csv")

# Save raw data

group_order <- c("CTRL","RZ","BZ","IZ","FZ")

walk2(mean_scores$gset, mean_scores$data, function(g,d) {
  
  b_plt <- d %>%
    ggplot(aes(x = factor(major_labl,
                          levels = group_order), y = score)) +
    geom_boxplot(outlier.size = 0) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          axis.text = element_text(size =  12),
          axis.title.x = element_text(size = 6),
          legend.position = "none") +
    ylab(g) +
    stat_compare_means(size = 4) +
    xlab("")
  
  pdf(paste0("./results/cell_death/wmean_score_", g, ".pdf"), height = 4, width = 4)
  plot(b_plt)
  dev.off()
  
})

pdf("./results/cell_death/n_nuclei.pdf", height = 4, width = 4)

b_plt <- pat_anns %>%
  ggplot(aes(x = factor(major_labl,
                        levels = group_order), y = ncells)) +
  geom_boxplot(outlier.size = 0) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size =  12),
        axis.title.x = element_text(size = 6),
        legend.position = "none") +
  ylab("No. nuclei") +
  stat_compare_means(size = 4,label.y = 15000) +
  xlab("")

plot(b_plt)

pat_anns %>% write_csv("./results/cell_death/n_nuclei.csv")

dev.off()
