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

liana_res <- readRDS("./results/cell_comms/FibMyeloid/liana_Fib_Myeloid.rds") %>% 
  liana_aggregate() %>%
  mutate(log10pvalue = -log10(cellphonedb.pvalue + 0.000001))

# Get the top 10 ligands 
myeloid <- c("SPP1_Macrophages", "LYVE_FOLR_Macrophages", 
             "LYVE_PLTP_Macrophages", "CCL18_Macrophages", 
             "Monocytes")

fibroblasts <- c("Fib_0", "Fib_3", "Myofib", "Fib_SCARA5")

# Myeloid interactions (senders)
myeloid_int <- liana_res %>%
  mutate(source_myeloid = source %in% myeloid,
         target_fibroblast = target %in% fibroblasts) %>%
  dplyr::filter(source_myeloid == TRUE,
                target_fibroblast == TRUE)

# For each combo we will get the top 10 interactions
myeloid_filt_ints <- myeloid_int %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:5) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

myeloid_int_plt <- myeloid_int %>%
  left_join(myeloid_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = myeloid_int_plt,
              source_groups = myeloid,
              target_groups = fibroblasts,
              specificity = 'log10pvalue')

pdf("./results/cell_comms/FibMyeloid/Myeloid_Fib.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# Fibroblast interactions (senders)
fib_int <- liana_res %>%
  mutate(source_fibroblast = source %in% fibroblasts,
         target_myeloid = target %in% myeloid) %>%
  dplyr::filter(source_fibroblast == TRUE,
                target_myeloid == TRUE)

# For each combo we will get the top 15 interactions
fib_filt_ints <- fib_int %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:5) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

fib_int_plt <- fib_int %>%
  left_join(fib_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(fib_int_plt,
                          source_groups = fibroblasts,
                          target_groups = myeloid,
                          specificity = 'log10pvalue')

pdf("./results/cell_comms/FibMyeloid/Fib_Myeloid.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# Focus now in the myeloid/fibroblast

# Myeloid interactions (senders)
myeloid_int <- liana_res %>%
  mutate(source_myeloid = source %in% myeloid,
         target_fibroblast = target %in% fibroblasts) %>%
  dplyr::filter(source_myeloid == TRUE,
                target_fibroblast == TRUE) %>%
  dplyr::filter(target == "Myofib")

# For each combo we will get the top 10 interactions
myeloid_filt_ints <- myeloid_int %>%
  dplyr::filter(source == "SPP1_Macrophages") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:35) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

myeloid_int_plt <- myeloid_int %>%
  left_join(myeloid_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(liana_agg = myeloid_int_plt,
                          source_groups = myeloid,
                          target_groups = fibroblasts,
                          specificity = 'log10pvalue')

pdf("./results/cell_comms/FibMyeloid/Myeloid_Myof.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()

# Fibroblast interactions (senders)
fib_int <- liana_res %>%
  mutate(source_fibroblast = source %in% fibroblasts,
         target_myeloid = target %in% myeloid) %>%
  dplyr::filter(source_fibroblast == TRUE,
                target_myeloid == TRUE) %>%
  dplyr::filter(target == "SPP1_Macrophages")

# For each combo we will get the top 15 interactions
fib_filt_ints <- fib_int %>%
  dplyr::filter(source == "Myofib") %>%
  dplyr::arrange(source, target, aggregate_rank) %>%
  group_by(source, target) %>%
  slice(1:35) %>%
  ungroup() %>%
  dplyr::select(ligand, receptor) %>%
  unique() %>%
  dplyr::mutate(keep = TRUE)

fib_int_plt <- fib_int %>%
  left_join(fib_filt_ints) %>%
  dplyr::filter(keep == TRUE)

test_plt <- liana_dotplot(fib_int_plt,
                          source_groups = fibroblasts,
                          target_groups = myeloid,
                          specificity = 'log10pvalue')

pdf("./results/cell_comms/FibMyeloid/SPP1_Fib.pdf", height = 6, width = 10)

plot(test_plt)

dev.off()
