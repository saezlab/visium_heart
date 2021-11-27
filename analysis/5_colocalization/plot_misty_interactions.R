# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generates cell maps that are generalizable

library(tidyverse)
library(cowplot)
library(mistyR)
library(Seurat)
source("./analysis/utils/misty_utilities.R")

path <- "./processed_visium/objects/"
outpath <- "./processed_visium/misty_red_objects/"
outpath_plots <- "./results/misty_red_objects_plts/"
slide_files <- list.files(path)
slide_names <- gsub("[.]rds", "", slide_files)
pat_ann <- readRDS("./markers/visium_patient_anns_revisions.rds")

param_df <- tibble(slide_name = slide_names,
                   visium_file = paste0(path, slide_files),
                   slide_out = paste0(outpath, slide_files %>% gsub("[.]rds", "_mistyassays.rds",.)),
                   plt_out = paste0(outpath_plots, slide_files %>% gsub("[.]rds", "_mistyassays.pdf",.)))

get_spat_contxt_plots <- function(slide_name, visium_file, slide_out, plt_out) {
  
  print(slide_name)
  
  visium_slide <- readRDS(visium_file)
  
  geometry <- GetTissueCoordinates(visium_slide,
                                   cols = c("row", "col"), 
                                   scale = NULL)
  
  data <- extract_seurat_data(visium.slide = visium_slide,
                              assay = "c2l",
                              geometry = geometry)
  
  #Juxta
  
  juxta_assay <- create_default_views(data = data,
                                      view.type = "juxta",
                                      view.param = 5,
                                      view.name = "juxta_c2l",
                                      spot.ids = rownames(data),
                                      geometry = geometry)
  
  juxta_assay <- juxta_assay$juxta_c2l_5$data %>%
    as.matrix()
  
  rownames(juxta_assay) <- rownames(data)
  
  visium_slide[["c2l_juxta"]] <- CreateAssayObject(data = t(juxta_assay))
  
  #Para
  
  para_assay <- create_default_views(data = data,
                                     view.type = "para",
                                     view.param = 15,
                                     view.name = "para_c2l",
                                     spot.ids = rownames(data),
                                     geometry = geometry)
  
  para_assay <- para_assay$para_c2l_15$data %>%
    as.matrix()
  
  rownames(para_assay) <- rownames(data)
  
  visium_slide[["c2l_para"]] <- CreateAssayObject(data = t(para_assay))
  
  # Plots
  
  cells <- rownames(GetAssayData(visium_slide, assay = "c2l"))
  
  DefaultAssay(visium_slide) <- "c2l"
  
  intra <- SpatialFeaturePlot(visium_slide, features = cells)
  
  DefaultAssay(visium_slide) <- "c2l_juxta"
  
  juxta <- SpatialFeaturePlot(visium_slide, features = cells)
  
  DefaultAssay(visium_slide) <- "c2l_para"
  
  para <- SpatialFeaturePlot(visium_slide, features = cells)
  
  DefaultAssay(visium_slide) <- "progeny"
  
  progeny <- SpatialFeaturePlot(visium_slide, features = rownames(visium_slide))
  
  pdf(file = plt_out, height = 15, width = 15)
  
  plot(intra)
  plot(juxta)
  plot(para)
  plot(progeny)
  
  dev.off()
  
  saveRDS(DietSeurat(
    visium_slide,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = c("c2l","c2l_juxta","c2l_para","progeny"),
    dimreducs = NULL
  ),
  slide_out)
  
  return(NULL)
  
}

pwalk(param_df, get_spat_contxt_plots)

# Create a guide to decide which slides to look at
misty_out_folder <- "./results/tissue_structure/misty/cell_map/"
misty_outs <- list.files(misty_out_folder, full.names = F)
misty_outs <- set_names(misty_outs, gsub("cm_", "", misty_outs) %>%
                          gsub("_c2l", "", .))

# Create a summary of MISTy results to navigate

misty_res <- collect_results(paste0(misty_out_folder, misty_outs))

R2_data <- misty_res$improvements %>%
  dplyr::filter(measure == "multi.R2") %>%
  dplyr::mutate(sample = gsub("_c2l", "", sample) %>%
                  strsplit(.,split = "cm_") %>%
                  map_chr(., ~ last(.x))) %>%
  dplyr::left_join(pat_ann, by = c("sample" = "sample_id")) %>%
  rename("R2" = value)

cell_order <- R2_data %>% 
  group_by(target) %>%
  summarize(med_value = median(R2)) %>%
  arrange(-med_value) %>%
  pull(target)

cells_R2_tile <- ggplot(R2_data, aes(x = factor(target,
                               levels = cell_order), 
                    y = sample, fill = R2)) +
 geom_tile() +
coord_equal() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cells_R2_box <- ggplot(R2_data, aes(x = factor(target,
                               levels = cell_order), y = R2)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
cells_R2_box_by_group <- ggplot(R2_data, aes(x = patient_group, y = R2)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(.~target)


write_csv(R2_data, "./results/tissue_structure/colocalization/c2l_R2_results.csv")

pdf("./results/tissue_structure/colocalization/c2l_R2_results.pdf", height = 7, width = 7)

plot(cells_R2_tile)
plot(cells_R2_box)
plot(cells_R2_box_by_group)

dev.off()

