# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Plot only fibroblasts and CM states vs niche
#' 
library(tidyverse)
library(ggpubr)
source("./analysis/utils/niche_utils.R")

pat_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")
ct_description <- readRDS("./results/ind_mats/cell_type_props.rds")
state_description_all <- readRDS("./results/ind_mats/cell_state_scorespos.rds")
spot_anns <- read_csv("./results/niche_mapping/mol_clust_class.csv")

int_res <- c("composition_niche",
             "Spatial_snn_res.0.2")

# CMs ---------------------------------------------------

res <- "Spatial_snn_res.0.2"
  
niche_info <- spot_anns %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
  dplyr::rename("mol_niche" = res) %>%
  dplyr::mutate(mol_niche = paste0("niche_", as.numeric(as.character(mol_niche)) + 1))

sel_niches <- c("niche_1", "niche_2", "niche_4")
states <- c("CM-damaged-CM", "CM-healthy-CM", "CM-intermediate-CM")

niche_info <- niche_info %>%
  dplyr::filter(mol_niche %in% sel_niches)

state_description <-  state_description_all %>%
  left_join(niche_info %>%
              dplyr::select(mol_niche, spot_id, patient_region_id)) %>%
  na.omit() %>%
  dplyr::filter(name %in% states)

my_comparisons <- list( c("niche_1", "niche_2"), c("niche_1", "niche_4"), c("niche_2", "niche_4") )

state_comparison_plt <- ggplot(state_description, aes(x = mol_niche, y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  facet_wrap(.~name, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized weighted mean") +
  xlab("")

pdf("./results/niche_mapping/Spatial_snn_res.0.2/CMSpatial_snn_res.0.2_bplot.pdf", height = 4, width = 6)

plot(state_comparison_plt)

dev.off()

# Fibs ---------------------------------------------------

res <- "composition_niche"
    
sel_niches <- c("niche_4", "niche_5")

niche_info <- spot_anns %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
  dplyr::rename("mol_niche" = res) %>%
  dplyr::mutate(mol_niche = paste0("niche_", mol_niche))

niche_info <- niche_info %>%
  dplyr::filter(mol_niche %in% sel_niches)

states <- state_description_all$name %>% unique()
states <- states[grepl("Fib", states) | grepl("Myeloid", states)]

state_description <-  state_description_all %>%
  left_join(niche_info %>%
              dplyr::select(mol_niche, spot_id, patient_region_id)) %>%
  na.omit() %>%
  dplyr::filter(name %in% states)

my_comparisons <- list( c("niche_4", "niche_5"))

state_comparison_plt <- ggplot(state_description, aes(x = mol_niche, y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  facet_wrap(.~name, scales = "free",nrow = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized weighted mean") +
  xlab("")

pdf("./results/niche_mapping/composition_niche/Fibcomposition_niche_bplot.pdf", height = 8, width = 10)

plot(state_comparison_plt)

dev.off()

# Endo

res <- "Spatial_snn_res.0.2"

sel_niches <- c("niche_10")

niche_info <- spot_anns %>%
  select_at(c("spot_id", "orig.ident", "patient_region_id", res)) %>%
  dplyr::rename("mol_niche" = res) %>%
  dplyr::mutate(mol_niche = paste0("niche_", as.numeric(as.character(mol_niche)) + 1))

niche_info <- niche_info %>%
  dplyr::filter(mol_niche %in% sel_niches)

states <- state_description_all$name %>% unique()
states <- states[grepl("Endo", states)]

state_description <-  state_description_all %>%
  left_join(niche_info %>%
              dplyr::select(mol_niche, spot_id, patient_region_id)) %>%
  na.omit() %>%
  dplyr::filter(name %in% states)

state_order <- state_description %>%
  group_by(name) %>%
  summarize(median_val= median(value)) %>%
  arrange(-median_val) %>%
  pull(name)

my_comparisons <- list(c("Endo-Capillary-Endo", "Endo-Venous-Endo"))

state_comparison_plt <- ggplot(state_description, aes(x = factor(name,
                                                                 levels = state_order), y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized weighted mean") +
  xlab("")

pdf("./results/niche_mapping/Spatial_snn_res.0.2/EndoSpatial_snn_res.0.2_bplot.pdf", height = 5, width = 4)

plot(state_comparison_plt)

dev.off()
