# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Differential gene expression between Cardiomyocytes 1 and 2 in CK160 (P3 BZ)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(ggrepel)


sample = 'CK160'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
DefaultAssay(sc) = 'RNA'


markers = FindMarkers(sc, ident.1 = 'Cardiomyocytes 1', ident.2 = 'Cardiomyocytes 2', 
	min.pct = 0.25, verbose = FALSE)
markers = markers[markers$p_val_adj < 0.05,]
markers$gene = rownames(markers)


write.table(markers[,c('gene',cols[-6])], 
	file = paste0(sample, '.CM1.vs.CM2.dge.txt'), quote = FALSE, sep = '\t', row.names = FALSE)



#---- Volcano plot of DGEs
# Courtesy of Zhijian Li (lzj1769)

df = read.table(file = paste0('../data/', sample, '.CM1.vs.CM2.dge.txt'), header = TRUE)
df$p_val_adj_log10 = -log10(df$p_val_adj) 
df$Sign = ifelse(-log10(df$p_val_adj) > 10 & abs(df$avg_logFC) > 0.3, 'Yes', 'No')


# Show interesting genes
df_text = subset(df, gene %in% c('ANKRD1', 'NPPB', 'MYO18B',
                                  'FHL2', 'CACNB2'))


pdf(file = 'cardiomyocyte_diff_genes_volcano.pdf', height = 6, width = 6)
ggplot(data = df, aes(x = avg_logFC, y = -log10(p_val_adj))) +
    geom_point(size = 1) +
    geom_label_repel(data = df_text, 
    	aes(x = avg_logFC, y = -log10(p_val_adj),
    		label = gene)) +
    xlab('log2 fold change') +
    ylab('-log10 p-value') +
    theme_cowplot()
dev.off()

