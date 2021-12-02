# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Here we try to identify doublet clusters from the state generator

library(SingleCellExperiment)
library(tidyverse)
library(HDF5Array)

ct_folder <- "./results/ct_data/"
cts <- list.dirs(ct_folder,recursive = F,full.names = F) %>%
  set_names()

map(cts, function(x) {
  print(x)
  sc_data_file <- paste0(ct_folder, x, "/", x, "_states_sce/")
  doublet_plt_file <- paste0(ct_folder, x, "/", x, "_doublets.pdf")
  sc_data <- loadHDF5SummarizedExperiment(sc_data_file)
  sc_meta <- colData(sc_data) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")
  pdf(doublet_plt_file, height = 4, width = 6)
  plot(ggplot(sc_meta, aes(x = opt_state, y = doublet_score)) +
    geom_boxplot())
  dev.off()
  
})



