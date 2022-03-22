# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Visualize trends of niche composition
library(tidyverse)
library(ComplexHeatmap)

pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

niche_props <- read_csv("results/niche_mapping/Spatial_snn_res.0.2/niche_props.csv") %>%
  left_join(pat_anns)

pw_test_res <- read_csv("results/niche_mapping/Spatial_snn_res.0.2/niche_props_pwtest_area_myogenic.csv")

niche_props <- niche_props %>%
  dplyr::filter(patient_group == "group_1",
                major_labl %in% c("BZ", "CTRL", "RZ"),
                mol_niche %in% pw_test_res$mol_niche) %>%
  group_by(mol_niche, major_labl) %>%
  summarize(mean_prop = mean(niche_prop)) %>%
  group_by(mol_niche) %>%
  dplyr::mutate(std_mean_prop = (mean_prop - mean(mean_prop))/sd(mean_prop)) %>%
  dplyr::filter(mol_niche %in% c("niche_0", "niche_1", "niche_3"))

myogenic_niche_plt <- niche_props %>%
  dplyr::select(-mean_prop) %>%
  pivot_wider(names_from = major_labl, values_from = std_mean_prop) %>%
  column_to_rownames("mol_niche") %>%
  as.matrix() %>%
  ComplexHeatmap::Heatmap(.,name = "std mean prop")

pdf("results/niche_mapping/Spatial_snn_res.0.2/niche_props_myogenic.pdf", height = 3,width = 4)

draw(myogenic_niche_plt)

dev.off()

niche_props %>% write_csv("results/niche_mapping/Spatial_snn_res.0.2/niche_props_myogenic.csv")
