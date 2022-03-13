# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Plotting imaging results

library(tidyverse)
library(ggpubr)

my_comparisons <- list( c("control", "MI"))

NPPB <- read_table2("./results/imaging/kuppe_rnascope.cm_states.count_table.2022_02_28.tsv") %>%
  ggplot(aes(x = group, y = NPPB_normalized)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))

ANKRD1 <- read_table2("./results/imaging/kuppe_rnascope.cm_states.count_table.2022_02_28.tsv") %>%
  ggplot(aes(x = group, y = ANKRD1_normalized)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=0.5))


pdf("./results/imaging/CM_NPPB_pos.pdf", height = 4, width = 2.5)

plot(NPPB)

dev.off()


pdf("./results/imaging/CM_ANKRD1_pos.pdf", height = 4, width = 2.5)

plot(ANKRD1)

dev.off()

