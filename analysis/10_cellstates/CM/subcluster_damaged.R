# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we cluster spots based on spot compositions
#' Aiming to explain stressed cardios

library(tidyverse)
library(cowplot)
library(Seurat)
library(factoextra)
library(compositions)
source("./analysis/utils/misty_pipeline.R")

path <- "./processed_visium/misty_red_objects/"
slide_files <- list.files(path)
slide_names <- gsub("_mistyassays[.]rds", "", slide_files)
outpath_plots <- "./results/stressed_CMs/"

param_df <- tibble(slide_name = slide_names,
                   visium_file = paste0(path, slide_files),
                   plt_out = paste0(outpath_plots, slide_files %>% gsub("[.]rds", ".pdf",.)))

# Here we identify clusters per slide
dissectCMs <- function(slide_name, visium_file, plt_out) {
  
  print(slide_name)
  
  visium_slide <- readRDS(visium_file) %>%
    positive_states(., assay = "cell_states") %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.10)
  
  stressed_scores <- GetAssayData(visium_slide, assay = "cell_states_pos")["CM-damaged-CM",]
  stressed_scores <- scale(stressed_scores[stressed_scores != 0])
  stressed_scores <- set_names(stressed_scores[,1], rownames(stressed_scores))
  
  #Let's take the top 5% of spots
  
  cut_off <- quantile(stressed_scores, 0.9)
  
  stressed_spots <- stressed_scores[stressed_scores >= cut_off] %>% names()
  
  # Let's filter the assay and cluster based on c2l
  
  visium_slide <- visium_slide[, stressed_spots]
  
  # Let's cluster now based on c2l
  
  c2l_mat <- GetAssayData(visium_slide, assay = "c2l_props") %>%
    t() 
  #%>%
   # scale()
  
  c2l_mat <- c2l_mat[, colnames(c2l_mat) %in%  c("Fib", "vSMCs", 
                                                 "Myeloid", "Endo", 
                                                 "Adipo")]
  
  # Generate ILR transformation
  baseILR <- ilrBase(x = c2l_mat,
                     method = "basic")
  cell_ilr <- as.matrix(ilr(c2l_mat, baseILR))
  colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))
  
  gex_hclust <- eclust(cell_ilr, "hclust", k = 3)
  
  clust_labels <- tibble(spot_id = names(gex_hclust$cluster),
                         label = paste0("clust_", gex_hclust$cluster)) %>%
    column_to_rownames("spot_id")
  
  visium_slide$opt_clust <- clust_labels[colnames(visium_slide),]
  
  # I will leverage the initial clustering information
  
  #Visualization
  
  DefaultAssay(visium_slide) <- "cell_states_pos"
  
  stressed_spatial <- SpatialFeaturePlot(visium_slide, features = "CM-damaged-CM", pt.size.factor = 2.)
  
  spatial_clust <- SpatialDimPlot(visium_slide, group.by = "opt_clust", pt.size.factor = 2)
  
  DefaultAssay(visium_slide) <- "c2l_props"
  
  c2l_spatial <- SpatialFeaturePlot(visium_slide, features = c("Fib", "vSMCs", 
                                                "Myeloid", "Endo", 
                                                "Adipo"), 
                     pt.size.factor = 1.8, 
                     combine = F)
  
  intra_plts <- cowplot::plot_grid(plotlist = c2l_spatial, ncol = 3, align = "hv")
  
  #DefaultAssay(visium_slide) <- "c2l_juxta"
  
  #c2l_spatial <- SpatialFeaturePlot(visium_slide, features = c("Fib", "vSMCs", 
  #                                                             "Myeloid", "Endo", 
  #                                                             "Adipo"), 
  #                                  pt.size.factor = 2.5, combine = F)
  
  #juxta_plts <- cowplot::plot_grid(plotlist = c2l_spatial, ncol = 3, align = "hv")
  
  c2l_vlns <- Seurat::VlnPlot(visium_slide,
                  features = c("Fib", "vSMCs", "Myeloid", "Endo", "Adipo"), 
                  assay = "c2l_props",
                  group.by = "opt_clust", combine = F)
  
  damaged_vln <- Seurat::VlnPlot(visium_slide,features = "CM-damaged-CM", assay = "cell_states_pos",
                  group.by = "opt_clust", combine = F)
  
  c2l_vlns <- cowplot::plot_grid(plotlist = c2l_vlns, ncol = 3, align = "hv")
  
  pdf(file= plt_out, width = 13, height = 10)
  
  plot(cowplot::plot_grid(stressed_spatial, spatial_clust, ncol = 2, align = "hv"))
  
  plot(intra_plts)
  
  #plot(juxta_plts)
  
  plot(cowplot::plot_grid(c2l_vlns, damaged_vln[[1]], ncol = 2, rel_widths = c(1, 0.5),align = "hv"))
  
  dev.off()
  
  return(list(mat = rbind(GetAssayData(visium_slide, assay = "c2l_props"),
               GetAssayData(visium_slide, assay = "cell_states_pos")["CM-damaged-CM",, drop = F]),
              annotations = visium_slide$opt_clust))
  
}

# Here analyze the individual results

stressed_res <- pmap(param_df, dissectCMs)

names(stressed_res) <- param_df$slide_name

# What about clustering spots from all slides? -------------------------------------------------------

all_spots <- map2(names(stressed_res), stressed_res, function(s_name, dat) {
  colnames(dat$mat) <- paste0(s_name, "..", colnames(dat$mat))
  return(dat$mat)
}) %>%
  reduce(., cbind) %>%
  t()

all_spots <- all_spots[, which(colnames(all_spots) %in% c("Fib", "vSMCs",
                                                          "Myeloid", "Endo",
                                                          "Adipo"))]

all_annotations <- map2(names(stressed_res), stressed_res, function(s_name, dat) {
  ann <- enframe(dat$annotations, name = "spot_id", value = "ind_label") %>%
    mutate(orig.ident = s_name)
  return(ann)
}) %>%
  enframe() %>%
  unnest() %>%
  dplyr::select(-name)

# Generate ILR transformation to do proper clustering

baseILR <- ilrBase(x = all_spots,
                   method = "basic")
cell_ilr <- as.matrix(ilr(all_spots, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

gex_hclust <- eclust(cell_ilr, "hclust", k = 3)

clust_labels <- tibble(spot_id = names(gex_hclust$cluster),
                       shared_label = paste0("sharedclust_", gex_hclust$cluster)) %>%
  mutate(orig.ident = strsplit(spot_id, "[.][.]") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(spot_id = strsplit(spot_id, "[.][.]") %>%
           map_chr(., ~.x[[2]]))

# Here we have the shared annotations

all_annotations <- left_join(all_annotations, clust_labels, by = c("spot_id", "orig.ident"))

# Here we merge proportions with annotations

props_df <- all_spots %>%
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  mutate(orig.ident = strsplit(spot_id, "[.][.]") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(spot_id = strsplit(spot_id, "[.][.]") %>%
           map_chr(., ~.x[[2]])) %>%
  left_join(all_annotations, by = c("spot_id", "orig.ident"))

# Proportions of shared 

# Annotate clusters with cell-types
# order = name (nest), clust_label, value

wilcox_wrap <- function(dat, clust_label = "ind_label") {
  niches <- dat[[clust_label]] %>%
    unique() %>%
    set_names()
  
  map(niches, function(g) {
    
    test_data <- dat %>%
      mutate(test_group = ifelse(.data[[clust_label]] == g,
                                 "target", "rest")) %>%
      mutate(test_group = factor(test_group,
                                 levels = c("target", "rest")))
    
    wilcox.test(value ~ test_group, 
                data = test_data,
                alternative = "greater") %>%
      broom::tidy()
  }) %>% enframe(clust_label) %>%
    unnest()
}

# Annotate shared clusters ---- 
# Panels and data

ct_annotation <- props_df %>%
  dplyr::select(-c("spot_id", "orig.ident")) %>%
  pivot_longer(-c("ind_label", "shared_label")) %>%
  dplyr::filter(name %in%  c("Fib", "vSMCs",
                             "Myeloid", "Endo",
                             "Adipo")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(wlcx_all = map(data, wilcox_wrap, clust_label = "shared_label")) %>%
  dplyr::select(wlcx_all) %>%
  unnest()

heatmap_wres <- ct_annotation %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value)) %>%
  dplyr::mutate(padj = ifelse(padj == 0, NA, padj)) %>%
  dplyr::mutate(padj = ifelse(is.na(padj), min(padj, na.rm = T), padj)) %>%
  dplyr::mutate(log_pvalue = -log10(padj)) %>%
  ggplot(., aes(x = name, y = shared_label, fill = log_pvalue)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
  ylab("") +
  xlab("") +
  coord_equal()

pdf("./results/stressed_CMs/sup_figures/cell_type_wilcox_tile.pdf", height = 4, width = 5)

plot(heatmap_wres)
ct_annotation %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value)) %>%
  dplyr::mutate(padj = ifelse(padj == 0, NA, padj)) %>%
  dplyr::mutate(padj = ifelse(is.na(padj), min(padj, na.rm = T), padj)) %>%
  dplyr::mutate(log_pvalue = -log10(padj)) %>%
  write_csv("./results/stressed_CMs/sup_figures/cell_type_wilcox_tile.csv")

dev.off()
  
bxplts_wres <- props_df %>%
  dplyr::select(-c("spot_id", "orig.ident")) %>%
  pivot_longer(-c("ind_label", "shared_label")) %>%
  dplyr::filter(name %in%  c("Fib", "vSMCs",
                             "Myeloid", "Endo",
                             "Adipo")) %>%
  ggplot(aes(x = shared_label, y = value)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  xlab("") +
  ylab("proportions") +
  facet_wrap(.~name, scales = "free_y") 

pdf("./results/stressed_CMs/sup_figures/cell_type_wilcox_box.pdf", height = 5, width = 5)

plot(bxplts_wres)

props_df %>%
  dplyr::select(-c("spot_id", "orig.ident")) %>%
  pivot_longer(-c("ind_label", "shared_label")) %>%
  dplyr::filter(name %in%  c("Fib", "vSMCs",
                             "Myeloid", "Endo",
                             "Adipo")) %>%
  write_csv("./results/stressed_CMs/sup_figures/cell_type_wilcox_box.csv")
  
dev.off()

# Here we identify:
# First the proportions of clusters per patient group
# Second, the slides to neglect in the last part

sample_anns <- read_csv("./markers/visium_patient_anns_revisions.csv")

nclusts_samples <- props_df %>%
  group_by(orig.ident, shared_label) %>%
  summarize(n_spots = n()) %>%
  dplyr::mutate(prop_clust = n_spots/sum(n_spots)) %>%
  left_join(sample_anns, by = c("orig.ident" = "sample_id"))
  

pdf("./results/stressed_CMs/sup_figures/pat_comp.pdf", height = 5, width = 5)

plot(nclusts_samples %>%
  ggplot(., aes(x = major_labl, y = prop_clust, color = major_labl)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  ggpubr::stat_compare_means(label.y = 1.1, size = 2) +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  facet_wrap(.~ shared_label, ncol = 4) +
  ylab("niche proportion"))

plot(nclusts_samples %>%
  ggplot(., aes(x = patient_group, y = prop_clust, color = patient_group)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  ggpubr::stat_compare_means(label.y = 1.1, size = 2) +
  facet_wrap(.~ shared_label, ncol = 4) +
  ylab("niche proportion"))

write_csv(nclusts_samples, "./results/stressed_CMs/sup_figures/pat_comp.csv")

dev.off()

# Check now which samples to drop from the other analysis

filtered_niches <- nclusts_samples %>% 
  summarise(n_clusts = n()) %>%
  dplyr::filter(n_clusts > 1) %>%
  pull(orig.ident)

# Perform kwallis of both clusters

kw_compilled <- map2(stressed_res[filtered_niches], filtered_niches, function(dat, s_name) {
  
  spc_annotations <- all_annotations %>%
    dplyr::filter(orig.ident == s_name)
  
  an_dat <- dat$mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spot_id") %>%
    left_join(spc_annotations, by = "spot_id") %>%
    dplyr::select(-c("spot_id", "orig.ident")) %>%
    pivot_longer(-c("ind_label", "shared_label")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(kw_res_ind = map(data, function(dt) {
      
      kruskal.test(value ~ ind_label, data = dt) %>%
        broom::tidy()
      
    })) %>%
    mutate(kw_res_all = map(data, function(dt) {
      
      if(length(unique(dt$shared_label)) > 0){
      
      kruskal.test(value ~ shared_label, data = dt) %>%
        broom::tidy()
        
      } else {
        NULL
      }
    })) %>%
    dplyr::select(kw_res_ind, kw_res_all)
  
})

kw_plt <- enframe(kw_compilled, name = "slide") %>% 
  unnest() %>%
  dplyr::select(-kw_res_ind) %>%
  unnest() %>% 
  dplyr::filter(name %in%  c("Fib", "vSMCs",
                             "Myeloid", "Endo",
                             "Adipo","CM-damaged-CM")) %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value)) %>%
  ggplot(aes(x = name, y = -log10(padj))) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  xlab("")

pdf( "./results/stressed_CMs/sup_figures/kw_results_perslide.pdf", height = 5, width = 5)

plot(kw_plt)

dev.off()

enframe(kw_compilled, name = "slide") %>% 
  unnest() %>%
  dplyr::select(-kw_res_ind) %>%
  unnest() %>% 
  dplyr::filter(name %in%  c("Fib", "vSMCs", 
                                          "Myeloid", "Endo", 
                                          "Adipo","CM-damaged-CM") ) %>%
  ungroup() %>%
  dplyr::mutate(padj = p.adjust(p.value)) %>%
  write_csv("./results/stressed_CMs/sup_figures/kw_results_perslide.csv")

# Repeat the plots with new annotation

# Here we identify clusters per slide
dissectCMs_v2 <- function(slide_name, visium_file, plt_out) {
  
  print(slide_name)
  
  visium_slide <- readRDS(visium_file) %>%
    positive_states(., assay = "cell_states") %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.10)
  
  stressed_scores <- GetAssayData(visium_slide, assay = "cell_states_pos")["CM-damaged-CM",]
  stressed_scores <- scale(stressed_scores[stressed_scores != 0])
  stressed_scores <- set_names(stressed_scores[,1], rownames(stressed_scores))
  
  #Let's take the top 5% of spots
  
  cut_off <- quantile(stressed_scores, 0.9)
  
  stressed_spots <- stressed_scores[stressed_scores >= cut_off] %>% names()
  
  # Let's filter the assay and cluster based on c2l
  
  visium_slide <- visium_slide[, stressed_spots]
  
  # Here you add new annotations
  clust_labels <- all_annotations %>%
    dplyr::filter(orig.ident == slide_name) %>%
    dplyr::select(-c("ind_label", "orig.ident")) %>%
    column_to_rownames("spot_id")
  
  visium_slide$opt_clust <- clust_labels[colnames(visium_slide), "shared_label"]
  
  # I will leverage the initial clustering information
  
  #Visualization
  
  DefaultAssay(visium_slide) <- "cell_states_pos"
  
  stressed_spatial <- SpatialFeaturePlot(visium_slide, features = "CM-damaged-CM", pt.size.factor = 2.)
  
  spatial_clust <- SpatialDimPlot(visium_slide, group.by = "opt_clust", pt.size.factor = 2)
  
  DefaultAssay(visium_slide) <- "c2l_props"
  
  c2l_spatial <- SpatialFeaturePlot(visium_slide, features = c("Fib", "vSMCs", 
                                                               "Myeloid", "Endo", 
                                                               "Adipo"), 
                                    pt.size.factor = 1.8, 
                                    combine = F)
  
  intra_plts <- cowplot::plot_grid(plotlist = c2l_spatial, ncol = 3, align = "hv")
  
  #DefaultAssay(visium_slide) <- "c2l_juxta"
  
  #c2l_spatial <- SpatialFeaturePlot(visium_slide, features = c("Fib", "vSMCs", 
  #                                                             "Myeloid", "Endo", 
  #                                                             "Adipo"), 
  #                                  pt.size.factor = 2.5, combine = F)
  
  #juxta_plts <- cowplot::plot_grid(plotlist = c2l_spatial, ncol = 3, align = "hv")
  
  c2l_vlns <- Seurat::VlnPlot(visium_slide,
                              features = c("Fib", "vSMCs", "Myeloid", "Endo", "Adipo"), 
                              assay = "c2l_props",
                              group.by = "opt_clust", combine = F)
  
  damaged_vln <- Seurat::VlnPlot(visium_slide,features = "CM-damaged-CM", assay = "cell_states_pos",
                                 group.by = "opt_clust", combine = F)
  
  c2l_vlns <- cowplot::plot_grid(plotlist = c2l_vlns, ncol = 3, align = "hv")
  
  pdf(file= plt_out, width = 13, height = 10)
  
  plot(cowplot::plot_grid(stressed_spatial, spatial_clust, ncol = 2, align = "hv"))
  
  plot(intra_plts)
  
  #plot(juxta_plts)
  
  plot(cowplot::plot_grid(c2l_vlns, damaged_vln[[1]], ncol = 2, rel_widths = c(1, 0.5),align = "hv"))
  
  dev.off()
  
  return(NULL)
  
}

stressed_res <- param_df %>%
  dplyr::mutate(plt_out = gsub("[.]pdf","_sharedclusts.pdf", plt_out)) %>%
  pmap(., dissectCMs_v2)

all_annotations %>%
  dplyr::select(-c("ind_label")) %>%
  write_csv("./results/stressed_CMs/sup_figures/spot_ann.csv")
