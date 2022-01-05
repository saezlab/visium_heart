# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we process results from the
#' Myeloid/Fibroblast interactions
#' 

library(Seurat)
library(liana)
library(tidyverse)

# Plot

liana_dotplot <- function(liana_agg,
                          source_groups,
                          target_groups,
                          specificity = "cellphonedb.pvalue",
                          magnitude = "sca.LRscore"){
  
  # Modify for the plot
  liana_mod <- liana_agg %>%
    # Filter to only the cells of interest
    filter(source %in% source_groups) %>%
    filter(target %in% target_groups) %>%
    rename(magnitude = !!magnitude) %>%
    rename(specificity = !!specificity) %>%
    unite(c("ligand", "receptor"), col = "interaction", sep = " -> ") %>%
    unite(c("source", "target"), col = "source_target", sep = " -> " ,remove = FALSE) %>%
    dplyr::select_at(c("source_target", "interaction", "magnitude", "specificity"))
  
  int_order <- liana_mod %>%
    arrange(source_target, -magnitude, -specificity) %>%
    pull(interaction) %>%
    unique()
  
  source_target_order <- liana_mod %>%
    arrange(source_target) %>%
    pull(source_target) %>%
    unique()
  
  liana_mod <- liana_mod %>%
    dplyr::mutate(interaction = factor(interaction, 
                                       levels = int_order),
                  source_target = factor(source_target,
                                         levels = source_target_order))
  
  # colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  
  
  # plot
  suppressWarnings(
    ggplot(liana_mod,
           aes(x = interaction,
               y = source_target,
               colour = magnitude,
               size = specificity
           )) +
      geom_point() +
      scale_color_gradientn(colours = viridis::viridis(20)) +
      theme_bw(base_size = 9) +
      theme(
        axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 9, colour = "gray6") #,
        # strip.text.y.left = element_text(angle = 0)
      ) +
      scale_y_discrete(position = "right") +
      labs(x = "Interactions (Ligand -> Receptor)",
           colour = "Expression\nMagnitude",
           size = "Interaction\nSpecificity",
           y = NULL
      ) +
      coord_equal()
  )
}

# Get results from myofibroblast and myeloid

liana_res <- readRDS("./results/cell_comms/Endo/liana_Endo.rds")

liana_res <- liana_res[which(names(liana_res) != "sca")] %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001))


# Get the top 10 ligands 
endos <- c( "Arterial_Endo", "Lymphatic_Endo", "Capillary_Endo", "Venous_Endo", "Endocardial_Endo")
cts <- c("Fib", "CM", "Endo", "vSMCs", "PC")

# PC (receivers)
state_reaceivers <- liana_res %>%
  mutate(source_cts= source %in% cts,
         target_states = target %in% endos) %>%
  dplyr::filter(source_cts == TRUE,
                target_states == TRUE)

# For each combo we will get the top 10 interactions
state_filt_ints <- state_reaceivers %>%
  dplyr::filter(source == "vSMCs",
                target == "Arterial_Endo") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:30) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_reaceivers_plt <- state_reaceivers %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = state_reaceivers_plt,
                          source_groups = "vSMCs",
                          target_groups = endos,
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/Endo/liana_arterialendo.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# Capillary
state_filt_ints <- state_reaceivers %>%
  dplyr::filter((source == "CM" | source == "PC"),
                target == "Capillary_Endo" ) %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:15) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_reaceivers_plt <- state_reaceivers %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = state_reaceivers_plt,
                          source_groups = c("CM", "PC"),
                          target_groups = endos,
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/Endo/liana_capillaryendo.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# PC (senders)
state_senders <- liana_res %>%
  mutate(source_states = source %in% endos,
         target_cts= target %in% cts) %>%
  dplyr::filter(source_states == TRUE,
                target_cts == TRUE)

# For each combo we will get the top 10 interactions
state_filt_ints <- state_senders %>%
  dplyr::filter(target == "vSMCs",
                source == "Arterial_Endo") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:30) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_senders_plt <- state_senders %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = state_senders_plt,
                          source_groups = endos,
                          target_groups = "vSMCs",
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/Endo/liana_arterialendo_sender.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# Capillary
# For each combo we will get the top 10 interactions
state_filt_ints <- state_senders %>%
  dplyr::filter((target == "CM" | target == "PC"),
                source == "Capillary_Endo") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:15) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_senders_plt <- state_senders %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = state_senders_plt,
                          source_groups = endos,
                          target_groups = c("CM", "PC"),
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/Endo/liana_capillaryendo_sender.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()
