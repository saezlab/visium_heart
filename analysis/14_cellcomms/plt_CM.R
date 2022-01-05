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

liana_res <- readRDS("./results/cell_comms/CM/liana_CM.rds")

cpdb <- liana_res$cellphonedb %>%
  dplyr::select(source, target, ligand, receptor, lr.mean)

liana_res <- liana_res[which(names(liana_res) != "sca")] %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001)) %>%
  dplyr::left_join(cpdb, by = c("source", "target", "ligand", "receptor"))

# Get the top 10 ligands 
cm <- c("damaged_CM")
cts <- c("Fib", "Adipo", "Endo", "vSMCs", "Myeloid")

# CM (receivers)
state_reaceivers <- liana_res %>%
  mutate(source_cts= source %in% cts,
         target_states = target %in% cm) %>%
  dplyr::filter(source_cts == TRUE,
                target_states == TRUE)

# For each combo we will get the top 10 interactions
state_filt_ints <- state_reaceivers %>%
  dplyr::filter(source %in% cts,
                target == "damaged_CM") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:10) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_reaceivers_plt_filtered <- state_reaceivers %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE) %>%
  group_by(ligand,receptor) %>%
  mutate(u_ints = n()) %>%
  ungroup() %>%
  dplyr::filter(u_ints < 3)

state_reaceivers_plt_unfiltered <- state_reaceivers %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE) %>%
  group_by(ligand,receptor) 

test_plt_rec_filtered <- liana_dotplot(liana_agg = state_reaceivers_plt_filtered,
                          source_groups = cts,
                          target_groups = cm,
                          specificity = 'log10pvalue',
                          magnitude = "lr.mean")

test_plt_rec_unfiltered <- liana_dotplot(liana_agg = state_reaceivers_plt_unfiltered,
                                       source_groups = cts,
                                       target_groups = cm,
                                       specificity = 'log10pvalue',
                                       magnitude = "lr.mean")

pdf("./results/cell_comms/CM/liana_SCMs_rec_filt.pdf", height = 6, width = 7)

plot(state_reaceivers_plt_filtered)

write_csv(state_reaceivers_plt_filtered, file = "./results/cell_comms/CM/liana_SCMs_rec_filt.csv")

dev.off()

pdf("./results/cell_comms/CM/liana_SCMs_rec_unfilt.pdf", height = 6, width = 10)

plot(test_plt_rec_unfiltered)

write_csv(state_reaceivers_plt_unfiltered, file = "./results/cell_comms/CM/liana_SCMs_rec_unfilt.csv")

dev.off()

# CM (senders)
state_senders <- liana_res %>%
  mutate(source_states = source %in% cm,
         target_cts= target %in% cts) %>%
  dplyr::filter(source_states == TRUE,
                target_cts == TRUE)

# For each combo we will get the top 10 interactions
state_filt_ints <- state_senders %>%
  dplyr::filter(target %in% cts,
                source == cm) %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:10) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

state_senders_plt_filt <- state_senders %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE) %>%
  group_by(ligand,receptor) %>%
  mutate(u_ints = n()) %>%
  ungroup() %>%
  dplyr::filter(u_ints < 3)

state_senders_plt_unfilt <- state_senders %>%
  left_join(state_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt_snd_filtered <- liana_dotplot(liana_agg = state_senders_plt_filt,
                          source_groups = cm,
                          target_groups = cts,
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

test_plt_snd_unfiltered <- liana_dotplot(liana_agg = state_senders_plt_unfilt,
                                       source_groups = cm,
                                       target_groups = cts,
                                       specificity = 'log10pvalue',
                                       magnitude = "logfc.logfc_comb")


pdf("./results/cell_comms/CM/liana_SCMs_snd_filt.pdf", height = 6, width = 7)

plot(test_plt_snd_filtered)

write_csv(state_senders_plt_filt, file = "./results/cell_comms/CM/liana_SCMs_snd_filt.csv")

dev.off()

pdf("./results/cell_comms/CM/liana_SCMs_snd_unfilt.pdf", height = 6, width = 10)

plot(test_plt_snd_unfiltered)

write_csv(state_senders_plt_unfilt, file = "./results/cell_comms/CM/liana_SCMs_snd_unfilt.csv")

dev.off()

