# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Miscellaneous feature/violin plots

library(Seurat)
library(ggplot2)
source('./sc_source/sc_source.R')


# Marker genes in CK158 (P1 C)
sample = 'CK158'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

genes = c('NPR3', 'CDH11', 'IL13RA')
for (gene in genes){
	pdf(file = paste0(sample, '_', gene, '.pdf'))
    print(VlnPlot(sc, features = gene, pt.size = 0, 
    	slot = 'scale.data', cols = cell.type.colors))
    dev.off()
} 



# Marker genes in integrated data set
sc = readRDS(file = '../data/scRNA.integrated.rds')
Idents(sc) = factor(Idents(sc), levels = cell.type.order)
DefaultAssay(sc)

genes = c('SEMA3G')
for (gene in genes){
	pdf(file = paste0('integrated_', gene, '.pdf'))
    print(VlnPlot(sc, features = gene, pt.size = 0, 
    	slot = 'scale.data', cols = cell.type.colors))
    dev.off()
} 

