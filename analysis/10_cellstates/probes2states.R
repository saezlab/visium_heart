# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we check to which states
#' the ordered probes belong

library(tidyverse)
library(cowplot)

probes <- read_csv("./markers/FibMyeloid_probe_genes.csv",col_names = F) %>%
  dplyr::rename("probe" = X1) %>%
  unique() %>%
  dplyr::pull(probe)

probes <- c(probes, "PLTP")

state_mrkrs <- tibble(marker_file = list.files("./cell_states", full.names = T)) %>%
  dplyr::mutate(cell_type = gsub("./cell_states/", "", marker_file)) %>%
  dplyr::mutate(marker_file = paste0(marker_file,"/annotation.rds")) %>%
  dplyr::mutate(markers = map(marker_file, readRDS)) %>%
  dplyr::select(cell_type, markers) %>%
  unnest() %>%
  dplyr::select(cell_type, p_val_adj, cluster, gene) %>%
  dplyr::filter(gene %in% probes,
                p_val_adj < 0.05) %>%
  dplyr::mutate(cluster = paste0(cell_type, "_", cluster)) %>%
  dplyr::mutate("ordered_probe" = T)


probe_plt <- ggplot(state_mrkrs, aes(y = factor(gene,
                                   levels = state_mrkrs$gene %>%
                                     unique()), x = cluster, fill = ordered_probe)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_equal() +
  ylab("") +
  xlab("")

pdf("./results/cell_states/probes2states.pdf", height = 7, width = 7)

plot(probe_plt)

dev.off()













