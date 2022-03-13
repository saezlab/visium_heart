# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generate plots of myeloid quantification
#' 

library(tidyverse)

# SPP1 is more abundant in IZ samples

spp1_postn <- read_table2("./results/spp1_quant_manual/kuppe_rnascope.fibroblast_myeloid.cell_marker_counts.2022_03_09.tsv")

my_comps <- list(c("control", "FZ"), c("control", "IZ"), c("FZ", "IZ"))

spp1_postn <- ggplot(spp1_postn, aes(x = sample_group, y = SPP1_CD163_norm_to_CD163, color = sample_group)) +
  geom_boxplot() +
  geom_point() +
  ylab("SPP1+ macrophage proportion") +
  ggpubr::stat_compare_means(comparisons = my_comps) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


pdf("./results/imaging/SPP1_quant_pos.pdf", width = 3,height = 4)

plot(spp1_postn)

dev.off()

# There are other cells besides SPP1

spp1_postn <- read_csv("./results/spp1_quant_manual/Mappe1_1.csv")

spp1_postn <- spp1_postn %>%
  mutate(sample = paste0("sample", seq(1, nrow(spp1_postn)))) %>%
  pivot_longer(-sample)


my_comps <- list(c("ccr2", "spp1"), c("ccr2", "trem2"), c("spp1", "trem2"))

trem2_plt <- ggplot(spp1_postn, aes(x = name, y = value)) +
  geom_boxplot() +
  geom_point(aes()) +
  ylab("Myeloid cell-type proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggpubr::stat_compare_means(paired = T,method = "wilcox.test", comparisons = my_comps)


pdf("./results/imaging/TREM2_quant_pos.pdf", width = 3,height = 4)

plot(trem2_plt)

dev.off()







