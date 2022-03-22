# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Functions associated to process liana results

library(tidyverse)

# Annotations in atlas
atlas_anns <- read_csv(file = "./processed_snrnaseq/cell_states/integrated_rnasamples_anns_wstates.csv")[,"annotation"] %>%
  pull() %>%
  unique()

# State markers
state_mrkrs <- tibble(marker_file = list.files("./cell_states", full.names = T)) %>%
  dplyr::mutate(cell_type = gsub("./cell_states/", "", marker_file)) %>%
  dplyr::mutate(marker_file = paste0(marker_file,"/annotation.rds")) %>%
  dplyr::mutate(markers = map(marker_file, readRDS)) %>%
  dplyr::select(cell_type, markers) %>%
  unnest() %>%
  dplyr::select(cell_type, p_val_adj, cluster, gene, avg_log2FC) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(-p_val_adj) %>%
  dplyr::filter(! grepl("Adipo", cell_type),
                ! grepl("Mast", cell_type)) %>%
  dplyr::select(cluster,gene) %>%
  dplyr::rename("cell_class" = cluster) %>%
  dplyr::mutate("right_ann" = cell_class %in% atlas_anns) %>%
  dplyr::mutate(cell_class = as.character(cell_class))

# Check wrong namings and manually change them

state_mrkrs %>% 
  dplyr::select(cell_class, right_ann) %>%
  unique() %>%
  dplyr::filter(right_ann == F)

state_mrkrs %>% 
  dplyr::select(cell_class, right_ann) %>%
  unique() %>%
  dplyr::filter(right_ann == T) %>%
  print(n = 50)

# Do the modifications needed

state_mrkrs <- state_mrkrs %>%
  mutate(cell_class = ifelse(cell_class == "DCs_FLT3_ITGAX", 
                             "Monocytes",
                             cell_class)) %>%
  mutate(cell_class = ifelse(cell_class == "Monocyte_CCL18", 
                             "CCL18_Macrophages",
                             cell_class)) %>%
  mutate(cell_class = ifelse(cell_class == "Monocyte_SPP1", 
                             "SPP1_Macrophages",
                             cell_class)) %>%
  dplyr::select(cell_class, gene)

# Plain cell-type markers

cell_type_mrkrs <- read_csv("./results/cell_markers/edgeR_cellmrkrs.csv") %>%
dplyr::filter(logFC > 0, FDR < 0.15) %>%
  arrange(name, FDR, - logFC) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  arrange(name, FDR, -logFC) %>%
  dplyr::select(name, gene) %>%
  dplyr::rename("cell_class" = name)

cell_class_mrkrs <- bind_rows(state_mrkrs, cell_type_mrkrs) %>%
  mutate("gene_exists" = TRUE)

  
# Function to plot LIANA results
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
                       filter_marker_genes = F,
                       out_dir) {
  
  print(source_groups)
  
  # (receivers)
  state_reaceivers <- liana_res %>%
    filter(source %in% source_groups,
           target %in% target_groups) %>%
    dplyr::filter(cellphonedb.pvalue <= (0.05))
  
  # Filter interactions based on markers
  
  if(filter_marker_genes == T) {
    
    state_reaceivers <- left_join(state_reaceivers, 
                      cell_class_mrkrs, 
                      by = c("source" = "cell_class", 
                             "ligand" = "gene")) %>%
      dplyr::rename("in_source" = gene_exists) %>%
      dplyr::filter(in_source == TRUE) %>%
      left_join(., 
                cell_class_mrkrs, 
                by = c("target" = "cell_class",
                       "receptor" = "gene")) %>%
      rename("in_target" = gene_exists) %>%
      dplyr::filter(in_target == TRUE)
      
  }
  
  # For each combo we will get the top 10 interactions
  # After filtering for specificity
  
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
                                         magnitude = "sca.LRscore")
  
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
