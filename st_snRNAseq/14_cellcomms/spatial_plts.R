# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate pairs of genes and pair of cells:

library(tidyverse)
library(Seurat)
library(viridis)
source("./analysis/utils/misty_pipeline.R")
source("./analysis/utils/spatial_plots_utils.R")

#' Need to first have a tibble containing the slide name
#' and pairs of cells
sample_dict <- read_csv("./markers/visium_patient_anns_revisions.csv")
slide_files_folder <- "./processed_visium/objects/"
slide_files <- list.files(slide_files_folder)
slide_ids <- gsub("[.]rds", "", slide_files)

slide_vector <- set_names(paste0(slide_files_folder, slide_files), slide_ids)

# We then make a function that loads each slide, filters by the
# genes defined, correlates them, returns a tibble

get_cor_df <-  function(f_path, gene_set, cts) {
  
  print(f_path)
  
  visium_slide <- readRDS(f_path)
  
  gene_set_filt <- gene_set[gene_set %in% rownames(visium_slide)]
  
  if(length(gene_set_filt) > 1) {
    
    useful_spots <- get_ct_spots(visium_slide, ct = cts, filter_prop = 0.1)
    visium_slide <- visium_slide[gene_set_filt, useful_spots]
    
    cor_markers <- cor(GetAssayData(visium_slide, assay = "SCT") %>%
                         as.matrix() %>%
                         t(), 
                       method = "spearman") %>%
      as.data.frame() %>%
      rownames_to_column("gene_a") %>%
      pivot_longer(-gene_a, names_to = "gene_b", values_to = "spearman_cor") %>%
      dplyr::filter(gene_a != gene_b) 
    
  } else {
    return(NULL)
  }
  
}

# Fib-Myeloid

interactions <- bind_rows(read_csv("./results/cell_comms/Fib/Myofib_rec__filt.csv"),
                      read_csv("./results/cell_comms/Fib/Myofib_snd__filt.csv")) %>%
  dplyr::select(source, ligand, target, receptor) %>%
  unique()

gene_set_list <- c(interactions$ligand, interactions$receptor) %>%
  unique()

fib_myel <- map(slide_vector, get_cor_df, gene_set = gene_set_list, cts = c("Fib", "Myeloid"))

# Get best in class:

#Myofib senders
send_bic <- fib_myel %>%
  enframe() %>%
  unnest() %>%
  left_join(interactions, by = c("gene_a" = "ligand", "gene_b" = "receptor")) %>%
  na.omit() %>%
  arrange(-spearman_cor) %>%
  group_by(gene_a, gene_b, source, target) %>%
  arrange(-spearman_cor) %>%
  dplyr::filter(spearman_cor > 0.15) %>%
  dplyr::filter(source == "Myofib",
                target == "SPP1_Macrophages")

#Myofib receivers
rec_bic <- fib_myel %>%
  enframe() %>%
  unnest() %>%
  left_join(interactions, by = c("gene_a" = "ligand", "gene_b" = "receptor")) %>%
  na.omit() %>%
  arrange(-spearman_cor) %>%
  group_by(gene_a, gene_b, source, target) %>%
  arrange(-spearman_cor) %>%
  dplyr::filter(spearman_cor > 0.15) %>%
  dplyr::filter(source == "SPP1_Macrophages",
                target == "Myofib")


plot_bics <- function(name, gene_a, gene_b, spearman_cor, source, target) {
  
  if(source == "SPP1_Macrophages") {
    
    source_alias = "Myeloid-Monocyte-SPP1"
    
  } else if (source == "Myofib") {
    
    source_alias = "Fib-Myofib"
    
  } 
  
  
  if(target == "SPP1_Macrophages") {
    
    target_alias = "Myeloid-Monocyte-SPP1"
    
  } else if (target == "Myofib") {
    
    target_alias = "Fib-Myofib"
    
  }
  
  visium_slide <- readRDS(paste0("./processed_visium/objects/", name, ".rds")) %>%
    positive_states(., assay = "cell_states") %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.1)
  
  
  DefaultAssay(visium_slide) <- "cell_states_pos"
  
  state_plot <- SpatialFeaturePlot(visium_slide, 
                                   features = c(source_alias, target_alias),
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 2) 
  
  DefaultAssay(visium_slide) <- "SCT"
  
  gene_plot <- SpatialFeaturePlot(visium_slide, 
                                   features = c(gene_a, gene_b),
                                   max.cutoff = "q99", 
                                   min.cutoff = "q1",
                                   stroke = 0,
                                   pt.size.factor = 2) 
  
  
  
  pdf(paste0("./results/state_structure/Fib_Myeloid/CCC_examples/",name,"_", source, "_", gene_a, "_", target, "_", gene_b, ".pdf"), height = 4, width = 8)
  
  plot(state_plot)
  plot(gene_plot)
  
  dev.off()
  
  
}

pwalk(rec_bic, plot_bics)

pwalk(send_bic, plot_bics)

