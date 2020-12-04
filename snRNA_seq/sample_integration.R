# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Integration of all samples with standard Seurat CCA integration pipeline

library(Seurat)
library(rlist)
library(pheatmap)
options(future.globals.maxSize = 30720*1024^2)
source('./sc_source/sc_source.R')

samples = c('CK158', 'CK159', 'CK160', 'CK162', 'CK163', 'CK164', 'CK165')
patients = c('P1_C', 'P5_CZ', 'P3_BZ', 'P4_CZ', 'P2_BZ',  'P3_RZ', 'P2_IZ')
new.id = as.data.frame(cbind(samples, patients))
rownames(new.id) = samples



obj.list = c()
for (sample in samples){
	sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
	sc$orig_annotation = Idents(sc)

	# Find variable features again
	sc = NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
	sc = FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
	sc = RenameCells(object = sc, add.cell.id = sample)
	obj.list = list.append(obj.list, sc)
}



# Integrate and cluster
res = 0.6
sc.anchors = FindIntegrationAnchors(object.list = obj.list, dims = 1:20, verbose = FALSE)
sc.int = IntegrateData(anchorset = sc.anchors, dims = 1:20, verbose = FALSE)
DefaultAssay(sc.int) = 'integrated'

sc.int = ScaleData(sc.int, verbose = FALSE)
sc.int = RunPCA(sc.int, npcs = 30, verbose = FALSE)
sc.int = RunUMAP(sc.int, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc.int = FindNeighbors(sc.int, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc.int = FindClusters(sc.int, resolution = res, verbose = FALSE)



# Add meta data 
sc.int$patient = new.id[sc.int$orig.ident, 'patients']
sc.int$patient_orig_annotation = paste(sc.int$patient, gsub(' ', '_', sc.int$orig_annotation), sep = '_')
Idents(sc.int) = sc.int$orig_annotation
Idents(sc.int) = factor(Idents(sc.int), levels = cell.type.order)



# Save data
saveRDS(sc.int, file = '../data/scRNA.integrated.rds')



#---- Visualise UMAP with cell annotations
pdf(file = paste0('scRNA_integrated_orig_annotation.pdf'), width = 12)
DimPlot(sc.int, reduction = 'umap', group.by = 'orig_annotation', cols = cell.type.colors)
dev.off()



#---- Subset and re-cluster fibroblast population for pseudotime analysis
sc = subset(sc.int, top_annotation == 'Fibroblasts')
DefaultAssay(sc) = 'integrated'

res = 0.5
sc = ScaleData(sc, verbose = FALSE)
sc = RunPCA(sc, npcs = 30, verbose = FALSE)
sc = RunUMAP(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc = FindNeighbors(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
sc = FindClusters(sc, resolution = res, verbose = FALSE)



# Save data 
saveRDS(sc, file = '../data/fib.subset.reclustered.rds')

