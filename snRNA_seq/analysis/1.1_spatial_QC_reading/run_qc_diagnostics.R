# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run the basic QC analysis of samples from their
#' spaceranger diagnostics
#' folder
#' |
#' sample---
#'          |
#'          ---spatial
#'          ---filtered_feature_bc_matrix.h5

library(tidyverse)
library(cowplot)

path <- "./visium_data/"

qc_feats <- c("sample_names",
              "Number of Spots Under Tissue",
              "Median Genes per Spot",
              "Mean Reads per Spot",
              "Fraction Reads in Spots Under Tissue",
              "Median UMI Counts per Spot",
              "Fraction of Spots Under Tissue")

sample_names <- list.files(path)

slide_files <- paste0(path,
                      sample_names,
                      "/outs/metrics_summary_csv.csv")

qc_stats <- tibble(sample_names,
                   qc_stats = map(slide_files, read_csv)) %>%
  unnest() %>%
  dplyr::select(all_of(qc_feats))

qc_stats_plts <- qc_stats %>%
  pivot_longer(-sample_names, names_to = "qc_feature") %>%
  group_by(qc_feature) %>%
  nest() %>%
  mutate(qc_plt = map2(qc_feature, data, function(dat_label, dat) {
    
    ggplot(dat, aes(x = sample_names,
                  y = value)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
            axis.text = element_text(size = 12),
            axis.title.y = element_text(size = 13))  +
      xlab("") + ylab(dat_label)
      
    
  }))

all_panels <- cowplot::plot_grid(plotlist = qc_stats_plts$qc_plt, align = "vh", ncol = 1)

pdf(height = 19, width = 20, file = "./processed_visium/initial_qc/all_qcs.pdf")

plot(all_panels)

dev.off()






















