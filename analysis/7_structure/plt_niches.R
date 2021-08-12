# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot "important" niches to space

library(tidyverse)
library(Seurat)
library(cowplot)

visium_folder = "./processed_visium/objects/"

sample_dict <- read.table("./markers/visium_annotations_ext.txt",
                          sep = "\t", header = T)

# Extract useful niches --------------------------------------

atlas_meta <- readRDS("./processed_visium/integration/ps_integrated_slides_niches.rds")[[1]][["annotations"]]

# Niche composition -----------------------------------
cell_info <- atlas_meta %>%
  dplyr::group_by(orig.ident, opt_clust_integrated) %>%
  summarize(ncells = length(opt_clust_integrated)) %>%
  dplyr::mutate(all_sample_cells = sum(ncells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cell_prop = ncells/all_sample_cells) %>%
  left_join(sample_dict, by = c("orig.ident" = "sample_id")) %>%
  mutate(niche = paste0("niche_", opt_clust_integrated))

# Niche composition in matrix
cell_props <- cell_info %>%
  dplyr::select(orig.ident, niche, cell_prop) %>%
  pivot_wider(values_from = cell_prop,
              names_from = niche) %>%
  column_to_rownames(var = "orig.ident") %>%
  as.data.frame()

cell_props[is.na(cell_props)] <- 0

head(cell_props)
dim(cell_props)

# Sort them as in the meta-dictionary
cell_props <- cell_props[sample_dict$sample_id,]
rownames(cell_props) <- sample_dict$pid

#Not all niches are informative for patient comparison(these are structures that are specific to a patient, noise, etc)
# So let's cut all the niches that aren't particularly present in 5 samples (minimum patient group size)
niche_filter <- colSums(cell_props > 0.01) >= 5
cell_props <- cell_props[, niche_filter]
niche_filter <- colnames(cell_props)

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample_name = visium_samples)

# Function that plots interesting clusters
plot_all_niches = function(visium_file, sample_name, cluster_list) {
  
  print(sample_name)
  
  visium_slide <- readRDS(visium_file)
  
  visium_slide <- subset(visium_slide, opt_clust_integrated %in% cluster_list)
  
  Idents(visium_slide) <- "opt_clust_integrated"

  niche_plt <-  SpatialDimPlot(visium_slide,
                 group.by = "opt_clust_integrated",
                 label = TRUE, 
                 label.size = 0,
                 stroke = 0, 
                 label.box = F) + ggtitle(sample_name)
  
  
  all_niches <- map(set_names(cluster_list), function(niche) {
    
    if(niche %in% Idents(visium_slide)) {
      SpatialDimPlot(visium_slide, 
                     cells.highlight = WhichCells(visium_slide, 
                                                  idents = niche),
                     label = F) + ggtitle(niche)
    } else {
      NULL
    }
  })
  
  all_niches <- plot_grid(plotlist = all_niches, 
                          ncol = 4, align = "hv",
                          labels = sample_name)
  
  
  pdf(file = paste0("./results/niche_mapping/mixed/niche_plts_", sample_name,".pdf"), width = 8, height = 9)
  
  plot(niche_plt)
  
  dev.off()
  
  pdf(file = paste0("./results/niche_mapping/sep/niche_plts_", sample_name,".pdf"), width = 18, height = 11)
  
  plot(all_niches)
  
  dev.off()
  
}


all_plts <- pwalk(visium_df, plot_all_niches, cluster_list = niche_filter)


# We load the slide meta-sata and map it














































