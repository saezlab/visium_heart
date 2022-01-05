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

liana_res <- readRDS("./results/cell_comms/PC/liana_PC.rds")

liana_res <- liana_res[which(names(liana_res) != "sca")] %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001))


# Get the top 10 ligands 
PCs <- c("Pericyte_1", "Pericyte_2")
cts <- c("Fib", "CM", "Endo", "vSMCs")

# PC (receivers)
PC_act <- liana_res %>%
  mutate(source_cts= source %in% cts,
         target_pcs = target %in% PCs) %>%
  dplyr::filter(source_cts == TRUE,
                target_pcs == TRUE)

# For each combo we will get the top 10 interactions
PC_filt_ints <- PC_act %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:10) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

PC_int_plt <- PC_act %>%
  left_join(PC_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = PC_int_plt,
                          source_groups = cts,
                          target_groups = "Pericyte_1",
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/PC/liana_PC.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# PC (senders)
PC_send <- liana_res %>%
  mutate(source_pcs= source %in% PCs,
         target_cts = target %in% cts) %>%
  dplyr::filter(source_pcs == TRUE,
                target_cts == TRUE)

# For each combo we will get the top 10 interactions
PC_filt_sends <- PC_send %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:10) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

PC_send_plt <- PC_send %>%
  left_join(PC_filt_sends) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = PC_send_plt,
                          source_groups = "Pericyte_1",
                          target_groups = cts,
                          specificity = 'log10pvalue',
                          magnitude = "logfc.logfc_comb")

pdf("./results/cell_comms/PC/liana_PC_senders.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()
