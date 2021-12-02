# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Estimate TF activities from 
#' a CM specific GRN

# Load libraries

library(tidyverse)
library(decoupleR)

# Load CM bulk contrast data

cm_data <- readRDS("./ext_data/tidy_de.RDS") %>%
  unnest() %>%
  dplyr::select(contrast, id, t) %>%
  pivot_wider(names_from = contrast, values_from = t) %>%
  column_to_rownames("id") %>%
  as.matrix()

# Load and filter GRN

regulons <- read_table2("./reg_nets/processed/CM.txt")

regulons <- regulons %>% 
  dplyr::filter(target %in% rownames(cm_data))

filtered_regulons <- regulons %>% 
  dplyr::select(source, target) %>%
  group_by(source) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]])) %>%
  mutate(n_targets = map_dbl(data, length)) %>%
  dplyr::filter(n_targets >= 5) %>%
  dplyr::select(-n_targets) %>%
  deframe() %>%
  names()

regulons <- regulons %>%
  dplyr::filter(source %in% filtered_regulons)

dorothea_res <- decoupleR::run_wmean(mat = cm_data, 
                                     network = regulons,
                                     .source = "source",
                                     .target = "target",
                                     .mor = "mor",
                                     .likelihood = "likelihood")

#

NFE2L1_res <- dorothea_res %>%
  dplyr::filter(source == "NFE2L1",
                statistic == "norm_wmean")
  
























