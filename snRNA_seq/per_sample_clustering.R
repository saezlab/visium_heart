# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Running standard Seurat clustering pipeline per sample individually

library(Seurat)
library(ggplot2)
source('./sc_source/sc_source.R')
options(future.globals.maxSize = 850*1024^2)
'%ni%' = Negate('%in%')
samples = c('CK158', 'CK159', 'CK160', 'CK161', 'CK162', 'CK163', 'CK164', 'CK165')
res = 0.5



for (sample in samples) {
	sc.data = Read10X(data.dir = paste0('../data/', sample, '/filtered_feature_bc_matrix'))
	sc = CreateSeuratObject(counts = sc.data, project = sample, min.cells = 3, min.features = 200)


	# Filter poor quality cells
	sc[['percent.mt']] = PercentageFeatureSet(sc, pattern = '^MT-')
	sc[['percent.ribo']] = PercentageFeatureSet(sc, pattern ="^RP[SL]")
	sc = subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & nCount_RNA > 300)


	# Cluster cells
	sc = NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
	sc = FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
	sc = ScaleData(sc, verbose = FALSE, features = rownames(sc))
	sc = RunPCA(sc, verbose = FALSE)
	sc = RunUMAP(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
	sc = FindNeighbors(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
	sc = FindClusters(sc, resolution = res, verbose = FALSE)


	# Save data
	saveRDS(sc, file = paste0('../data/', sample, '.rds'))


	# Get gene specificity score and conditionial probability of expression per cluster
	# Use for cluster annotations
	sg = run_genesorter(sc, write.file = TRUE, file.name.prefix = sample)


	# Cluster quality plots
	pdf(file = paste0(sample, '_cluster_quality.pdf'), width = 10)
	features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo')
	for (feature in features){
		print(FeaturePlot(sc, features = feature, reduction = 'umap') + 
			labs(title = sample, colour = feature) +  
			theme(legend.title = element_text(size = 12)))
	}
	dev.off()
}



#---- Filter poor quality clusters and add annotations

# CK158 (P1 C)
sample = 'CK158'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('8', '12', '13', '15', '16', '17', '18', '19')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Fibroblasts 1', `1` = 'Macrophages 1',
    `2` = 'Endothelial cells 1', `3` = 'Cardiomyocytes', `4` = 'Fibroblasts 2', 
    `5` = 'Pericytes', `6` = 'Endothelial cells 2', `7` = 'Endothelial cells 3',
    `9` =  'Macrophages 2', `10` = 'T cells', `11` = 'Neuronal cells',
    `14` = 'Vascular smooth muscle cells')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK159 (P5 CZ)
sample = 'CK159'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('3', '5', '11', '13', '14')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Fibroblasts 1', `1` = 'Endothelial cells 1',
    `2` = 'Pericytes', `4` = 'Macrophages 1', `6` = 'Cardiomyocytes', 
    `7` = 'Endothelial cells 2', `8` = 'Endothelial cells 4',
    `9` = 'Endothelial cells 3', `10` = 'Neuronal cells', `12` = 'T cells')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK160 (P3 BZ)
sample = 'CK160'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('1', '12')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Fibroblasts 1', `2` = 'Cardiomyocytes 1',
    `3` = 'Macrophages 1', `4` = 'Endothelial cells 3', `5` = 'Fibroblasts 2',
    `6` = 'Cardiomyocytes 2', `7` = 'Pericytes', `8` = 'T cells',
    `9` = 'Endothelial cells 1', `10` = 'Neuronal cells',
    `11` = 'Vascular smooth muscle cells', `13` = 'Lymphatic endothelial cells')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK161 (P3 IZ) 
# Discard entire sample, too poor quality



# CK162 (P4 CZ)
sample = 'CK162'
sc = readRDS(file = paste0('../data/', sample, '.rds'))

sc = RenameIdents(sc, `0` = 'Endothelial cells 1', `1` = 'Fibroblasts 3', 
    `2`= 'Endothelial cells (damaged)', `3` = 'Endothelial cells 3',
    `4` = 'Endothelial cells 4', `5` = 'Endothelial cells 5', `6` = 'Pericytes',
    `7` = 'Vascular smooth muscle cells', `8` = 'Fibroblasts 4',
    `9` = 'Cardiomyocytes', `10` = 'Fibroblasts 5', `11` = 'Macrophages 1',
    `12` = 'T cells', `13` = 'Neuronal cells', `14` = 'Lymphatic endothelial cells')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK163 (P2 BZ)
sample = 'CK163'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('7', '12', '13', '15', '16')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Fibroblasts 1', `1` = 'Endothelial cells 1',
    `2` = 'Pericytes', `3` = 'Cardiomyocytes', `4` = 'Macrophages 1',
    `5` = 'Endothelial cells 1', `6` = 'Endothelial cells 3', 
    `8` = 'Fibroblasts 2', `9` = 'Vascular smooth muscle cells', 
    `10` = 'T cells', `11` = 'Neuronal cells', `14` = 'Erythrocytes')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK164 (P3 RZ)
sample = 'CK164'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('4', '12', '13')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `0` = 'Fibroblasts 1', `1` = 'Cardiomyocytes',
    `2` = 'Macrophages 1', `3` = 'Fibroblasts 2', `5` = 'Endothelial cells 1', 
    `6` = 'Pericytes', `7` = 'T cells', `8` = 'Endothelial cells 2', 
    `9` = 'Lymphatic endothelial cells', `10` = 'Neuronal cells', `11` = 'Adipocytes')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



# CK165 (P2 IZ)
sample = 'CK165'
sc = readRDS(file = paste0('../data/', sample, '.rds'))
bad.clusters = c('0', '4', '6')
sc = subset(sc, subset = seurat_clusters %ni% bad.clusters)

sc = RenameIdents(sc, `1` = 'Cardiomyocytes', `2` = 'Endothelial cells 1', 
    `3` = 'Fibroblasts 0', `5` = 'Pericytes')
sc$top_annotation = gsub(' [0-9]', '', Idents(sc))
saveRDS(sc, file = paste0('../data/', sample, '.filtered.annotated.rds'))



#---- Visualise UMAPs with cell annotations
for (sample in samples){
	sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
	Idents(sc) = factor(Idents(sc), levels = cell.type.order)


	pdf(file = paste0(sample, '.pdf'), height = 5)
	print(DimPlot(sc, reduction = 'umap', cols = cell.type.colors))
	dev.off()
}

