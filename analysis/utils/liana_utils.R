# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Functions associated to process liana results

library(tidyverse)

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

# This function automatically filters 
liana_outs <- function(source_groups,
                       target_groups,
                       top = 10,
                       file_alias,
                       max_interactions = 2,
                       out_dir) {
  
  print(source_groups)
  
  # Endo (receivers)
  state_reaceivers <- liana_res %>%
    filter(source %in% source_groups,
           target %in% target_groups)
  
  # For each combo we will get the top 10 interactions
  state_filt_ints <- state_reaceivers %>%
    dplyr::arrange(source, target, aggregate_rank) %>%
    group_by(source, target) %>%
    slice(1:top) %>%
    ungroup() %>%
    dplyr::select(ligand, receptor) %>%
    unique() %>%
    dplyr::mutate(keep = TRUE)
  
  # We map to the original 
  state_reaceivers_plt_filtered <- state_reaceivers %>%
    left_join(state_filt_ints) %>%
    dplyr::filter(keep == TRUE) %>%
    dplyr::mutate(specific = log10pvalue > 3) %>%
    group_by(ligand,receptor) %>%
    mutate(u_ints = sum(specific)) %>%
    ungroup() %>%
    dplyr::filter(u_ints <= max_interactions)
  
  state_reaceivers_plt_unfiltered <- state_reaceivers %>%
    left_join(state_filt_ints) %>%
    dplyr::filter(keep == TRUE) %>%
    group_by(ligand,receptor) 
  
  test_plt_rec_filtered <- liana_dotplot(liana_agg = state_reaceivers_plt_filtered,
                                         source_groups = source_groups,
                                         target_groups = target_groups,
                                         specificity = 'log10pvalue',
                                         magnitude = "lr.mean")
  
  pdf(paste0(out_dir, file_alias,"_filt.pdf"), height = 6, width = 7)
  
  plot(test_plt_rec_filtered)
  
  write_csv(state_reaceivers_plt_filtered, 
            file = paste0(out_dir, file_alias,"_filt.csv"))
  
  dev.off()
  
  test_plt_rec_unfiltered <- liana_dotplot(liana_agg = state_reaceivers_plt_unfiltered,
                                           source_groups = source_groups,
                                           target_groups = target_groups,
                                           specificity = 'log10pvalue',
                                           magnitude = "lr.mean")
  
  pdf(paste0(out_dir, file_alias,"_unfilt.pdf"), height = 6, width = 10)
  
  plot(test_plt_rec_unfiltered)
  
  write_csv(state_reaceivers_plt_unfiltered, 
            file = paste0(out_dir, file_alias,"_unfilt.csv"))
  
  dev.off()
  
}
