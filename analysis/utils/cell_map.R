# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we compare the markers of niches coming
#' from similar structures

library(Seurat)
library(tidyverse)

# Here per spot we have the cell-tye
ct_description <- readRDS("./results/ind_mats/cell_type_props.rds")

# Let's keep maximum proportion
ct_description <- ct_description %>%
  group_by(spot_id) %>%
  dplyr::filter(value == max(value))

# If you aren't Fib/CM/Myeloid, you are other
plt_cts <- c("CM", "Fib", "Myeloid")

ct_description <- ct_description %>%
  dplyr::mutate(name = ifelse(name %in% plt_cts, name, "other")) %>%
  dplyr::mutate(orig.ident = strsplit(spot_id, "[.][.]") %>%
                  map_chr(., ~ .x[1]),
                raw_id = strsplit(spot_id, "[.][.]") %>%
                  map_chr(., ~ .x[2])) %>%
  dplyr::select(-c("spot_id", "value")) %>%
  dplyr::rename("plt_lbl" = name)

# Here we plot the maximum

visium_folder = "./processed_visium/objects/"

colors <- c("#c33124", "#a1dffb", "#e8a628", "darkgrey")

# Get visium slides --------------------------------
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample_name = visium_samples)

# Function that plots interesting clusters
plot_all_niches = function(visium_file, sample_name) {
  
  print(sample_name)
  
  visium_slide <- readRDS(visium_file)
  
  meta_data <- visium_slide@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("raw_id") %>%
    left_join(ct_description, by = c("orig.ident", "raw_id"))
  
  visium_slide$plt_lbl <- meta_data$plt_lbl
  
  niche_plt <-  SpatialDimPlot(visium_slide,
                               group.by = "plt_lbl",
                               label = TRUE, 
                               label.size = 0,
                               stroke = NA, 
                               label.box = F,) + 
    scale_fill_manual(values = colors)
  
  pdf(paste0("./results/cell_type_map/", sample_name, ".pdf"), height = 3, width = 3)
  
  print(niche_plt)
  
  dev.off()

}

pwalk(visium_df, plot_all_niches)



group_by(ct_description)













