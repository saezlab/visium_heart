# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we rgenerate basic QC plots of an integrated object
#' These include:
#' -Number of cells
#' -Distribution of counts
#' -Distribution of features
#' -Distribution of percentage of mitochondrial genes

library(Seurat)
library(tidyverse)
library(cowplot)
library(optparse)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "where's your Seurat object?"),
  make_option(c("--out_fig_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the fig objects"))

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Main ------------------------------------------------------------------

scell <- readRDS(path)

meta_data <- scell@meta.data

ncells_p <- meta_data %>%
  group_by(orig.ident) %>%
  summarise(ncells = length(orig.ident)) %>%
  ggplot(aes(x = orig.ident, y = ncells)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of spots")

nCount_RNA_p <- meta_data  %>%
  ggplot(aes(x = orig.ident, y = nCount_Spatial)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of UMIs")

nFeature_RNA_p <-  meta_data  %>%
  ggplot(aes(x = orig.ident, y = nFeature_Spatial)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Number of genes")

percent_mt_p <- meta_data  %>%
  ggplot(aes(x = orig.ident, y = percent.mt)) +
  geom_violin() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank()) +
  ylab("Mitochondrial gene percentage")

# Generate panel

panel_p <- cowplot::plot_grid(ncells_p, nCount_RNA_p, 
                              nFeature_RNA_p, percent_mt_p,
                              nrow = 2, ncol = 2, align = "hv")

umap_ps <- FeaturePlot(scell, features = c("nCount_Spatial", 
                                             "nFeature_Spatial",
                                             "percent.mt"),
                   combine = F, raster = T)

pdf(out_fig_path, height = 7, width = 7)

print(panel_p)

walk(umap_ps, print)

dev.off()

