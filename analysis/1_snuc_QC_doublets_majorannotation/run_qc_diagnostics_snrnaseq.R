# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run the basic QC analysis of samples from their
#' spaceranger diagnostics
#' folder
#' |
#' sample---outs
#'          |
#'          ---spatial
#'          ---filtered_feature_bc_matrix.h5

library(tidyverse)
library(cowplot)

path <- "./snrnaseq_data/"

qc_feats <- c("sample_names",
              "Estimated Number of Cells",
              "Mean Reads per Cell",
              "Median Genes per Cell",
              "Fraction Reads in Cells",
              "Total Genes Detected",
              "Median UMI Counts per Cell")

sample_names <- list.files(path)

slide_files <- paste0(path,
                      sample_names,
                      "/outs/metrics_summary.csv")

qc_stats <- tibble(sample_names,
                   qc_stats = map(slide_files, read_csv)) %>%
  unnest() %>%
  dplyr::select(all_of(qc_feats))

qc_stats$`Fraction Reads in Cells` <- gsub("%", "", qc_stats$`Fraction Reads in Cells`) %>%
  as.numeric()

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

pdf(height = 20, width = 17, file = "./processed_snrnaseq/initial_qc/all_qcs.pdf")

plot(all_panels)

dev.off()
