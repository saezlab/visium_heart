library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 30720*1024^2)

celltypes_2_top <- c("Cardiomyocytes 1" = "Cardiomyocytes",
                     "Cardiomyocytes 2" = "Cardiomyocytes",
                     "Endothelial cells 1" = "Endothelial cells",
                     "Endothelial cells 2" = "Endothelial cells",
                     "Endothelial cells 3" = "Endothelial cells",
                     "Endothelial cells 4" = "Endothelial cells",
                     "Fibroblasts 0" = "Fibroblasts",
                     "Fibroblasts 1" = "Fibroblasts",
                     "Fibroblasts 2" = "Fibroblasts",
                     "Fibroblasts 3" = "Fibroblasts",
                     "Fibroblasts 4" = "Fibroblasts",
                     "Fibroblasts 5" = "Fibroblasts",
                     "Macrophages 1" = "Macrophages",
                     "Macrophages 2" = "Macrophages")

samples_2_patients <- c("CK166" = "P1_RZ",
                        "CK167" = "P5_CR",
                        "CK168" = "P3_BZ",
                        "CK170" = "P4_CR",
                        "CK171" = "P2_BZ",
                        "CK173" = "P3_RZ",
                        "CK174" = "P2_IZ")

cells.col <- c("Cardiomyocytes" = "#800000",
               "Fibroblasts" = "#000075",
               "Endothelial cells" = "#e6194B",
               "Pericytes" = "#f58231",
               "T cells" = "#ffe119",
               "Lymphatic endothelial cells" = "#ffd8b1",
               "Adipocytes" = "#000000",
               "Neuronal cells" = "#42d4f4",
               "Erythrocytes" = "#D3D3D3",
               "Damaged endothelial cells" = "#999999",
               "Vascular smooth muscle cells" = "#aaffc3")

obj <- readRDS("../data/heart.integrated.Rds")
df <- as.data.frame(obj@meta.data)

df$top_annotation <- stringr::str_replace_all(df$celltype, celltypes_2_top)
df$patient <- stringr::str_replace_all(df$patient.ident, samples_2_patients)

p <- ggplot(data = df, aes(x = patient, fill = as.factor(top_annotation))) +
    geom_bar(position = position_fill(reverse = TRUE)) + 
    labs(y = 'Cell type proportion', x = element_blank(), fill = 'Cell type') + 
    theme_classic() + scale_fill_manual(values = my.colors)


pdf(file = '../data/cell_proportions.pdf')
print(p)
dev.off()