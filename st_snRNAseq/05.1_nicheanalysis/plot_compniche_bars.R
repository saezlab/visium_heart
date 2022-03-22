# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# In this script I will generate the stacked barplots of the compositional niches


library(Seurat)
library(tidyverse)
library(viridis)
source("./analysis/utils/spatial_plots_utils.R")


niche_cols <- list(niche_1 = "#D51F26",
                   niche_2 = "#272E6A",
                   niche_3 = "#208A42",
                   niche_4 = "#89288F",
                   niche_5 = "#F47D2B",
                   niche_6 = "#FEE500",
                   niche_7 = "#8A9FD1",
                   niche_8 = "#C06CAB",
                   niche_9 = "#D8A767") %>%
  unlist()


niche_props <- "./results/niche_mapping/composition_niche/niche_props.csv" %>%
  read_csv()

niche_props_plt <- ggplot(niche_props,aes(x = patient_region_id, y = niche_prop, fill = mol_niche)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 10, 
                                   hjust = 1),
        axis.text.y = element_text(size = 10)) +
  scale_fill_manual(values = niche_cols)  +
  ylab("") +
  xlab("")

pdf("./results/niche_mapping/composition_niche/niche_props_bars.pdf", height = 4, width = 7)

plot(niche_props_plt)

dev.off()












