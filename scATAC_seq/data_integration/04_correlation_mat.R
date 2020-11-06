library(Seurat)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

obj <- readRDS("../data/heart.integrated.Rds")

Idents(obj) <- "patient_celltype"
df <- Reduce("rbind", AverageExpression(object = obj,
                                        assays="integrated"))

celltype_cols = c('Cardiomyocytes' = '#800000',
                  'Cardiomyocytes_1' = '#800000',
                  'Cardiomyocytes_2' = '#9A6324',
                  'Fibroblasts_1' = '#911eb4',
                  'Fibroblasts_2' = '#e6beff',
                  'Fibroblasts_3' = '#f032e6',
                  'Fibroblasts_4' = '#f032e6',
                  'Fibroblasts_5' = '#f032e6',
                  'Fibroblasts_0' = '#f032e6',
                  'Endothelial_cells_1' = '#000075',
                  'Endothelial_cells_2' = 'blue',
                  'Endothelial_cells_3' = '#568198',
                  'Endothelial_cells_4' = '#469990',
                  'Endothelial_cells_5' = '#469990',
                  'Endothelial_cells_6' = '#469990',
                  'Macrophages' = '#e6194B',
                  'Macrophages_1' = '#e6194B',
                  'Macrophages_2' = '#fabebe',
                  'Pericytes' = '#f58231',
                  'T_cells' = '#ffe119',
                  'Lymphatic_endothelial_cells' = '#ffd8b1',
                  'Adipocytes' = '#000000',
                  'Neuronal_cells' = '#42d4f4',
                  'Erythrocytes' = '#fdbf6f',
                  'Endothelial_cells_(damaged)' = '#999999',
                  'Vascular_smooth_muscle_cells' = '#aaffc3')

df.plot <- cor(df)

patients <- stringr::str_split_fixed(colnames(df), "_", 3)[, 1]
locations <- stringr::str_split_fixed(colnames(df), "_", 3)[, 2]
celltypes <- stringr::str_split_fixed(colnames(df), "_", 3)[, 3]

celltypes <- stringr::str_replace_all(group,  " ", "_")

col <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

ha = HeatmapAnnotation(
    patients = patients,
    locations = locations,
    celltypes = celltypes, 
    col = list(patients = c("P1" = "#e41a1c",
                            "P2" = "#377eb8",
                            "P3" = "#4daf4a",
                            "P4" = "#984ea3",
                            "P5" = "#ff7f00"),
               locations = c("RZ" = "#ffff33",
                             "CR" = "#a65628",
                             "BZ" = "#f781bf",
                             "IZ" = "#b2df8a"),
               celltypes = celltype_cols
    ),
    show_annotation_name = FALSE
)

p <- Heatmap(matrix = df.plot,
             col = col,
             name = "Pearson",
             show_row_names = FALSE,
             show_column_names = FALSE,
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "pearson",
             clustering_method_columns = "average",
             clustering_method_rows = "average",
             cluster_rows = cluster_within_group(df.plot, celltypes),
             cluster_columns = cluster_within_group(df.plot, celltypes),
             top_annotation = ha)

pdf("scATAC_integrated_patient_celltype_correlation_v2.pdf", width = 10, height = 8)
draw(p)
dev.off()
