# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Pseudotime analysis with Monocle 3 on integrated fibroblasts

library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(gprofiler2)
library(gplots)
library(mclust)
library(genesorteR)
library(monocle3)
library(dplyr)
source('./sc_source/sc_source.R')


# Pseudotime analysis with Monocle3
# Convert dataset to monocle3 cell dataset object
sc = readRDS(file = '../data/fib.subset.reclustered.rds')
cds = as.cell_data_set(sc)


# Cluster cells with Leiden clustering
cds = cluster_cells(cds, cluster_methods = 'leiden', reduction_method = 'UMAP')


sc.sub = subset(as.Seurat(cds), monocle3_partitions == 1)
cds = as.cell_data_set(sc.sub)


# Trajectory analysis
cds = learn_graph(cds)


# Get root cells in trajectory, here most SCARA5+
# Order cells along pseudotime
root.gene = 'SCARA5'
max.root = which.max(unlist(FetchData(sc.sub, root.gene)))
max.root = colnames(sc.sub)[max.root]
cds = order_cells(cds, root_cells = max.root)


sc.sub = as.Seurat(cds, assay = 'integrated')
sc$monocle3_clusters = sc.sub$monocle3_clusters
sc$monocle3_partitions = sc.sub$monocle3_partitions
sc$monocle3_pseudotime = sc.sub$monocle3_pseudotime


pdf(file = 'fib_monocle3_pseudotime.pdf')
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
	label_branch_points = FALSE)
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, 
	label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, 
	label_leaves = FALSE, label_branch_points = FALSE,
	show_trajectory_graph = FALSE)
FeaturePlot(sc.sub, feature = 'monocle3_pseudotime')
dev.off()



#---- Differentially expressed genes along the trajectory

cds = estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[['gene_short_name']] = rownames(sc[['RNA']])
cds_dif_test = graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

write.table(cds_dif_test[,c(5,1,2,3,4,6)], file = 'fib.monocle3.dge.genes.txt', 
	sep = '\t', quote = FALSE, row.names = FALSE)



# Visualise top 40 genes in heatmap ordered on pseudotime
n = 40
cds_dif_test = cds_dif_test[cds_dif_test$q_value < 0.05,]

gene.list = cds_dif_test %>%
	arrange(desc(morans_test_statistic), desc(-q_value)) %>% 
	rownames()

top.genes = gene.list[1:n]


# Cluster cells (based on pseudotime) for heatmap
kk = Mclust(-sc$monocle3_pseudotime, 5, modelNames = 'E')


# Cluster genes with k-means clustering
plot = plotMarkerHeat(sc@assays$RNA@data[,order(sc$monocle3_pseudotime, decreasing = FALSE)], 
	kk$classification[order(sc$monocle3_pseudotime, decreasing = TRUE)], 
	top.genes, averageCells = 10^1, clusterGenes = TRUE, clusterGenesK = 3, 
	gap = FALSE, outs = TRUE, plotheat = FALSE)


# Change ordering of gene clusters manually so SCARA5+ top, POSTN+ bottom 
gene_order = plot$gene_class_info
gene_order[plot$gene_class_info == 3] = 2
gene_order[plot$gene_class_info == 2] = 3


pdf(file = 'fib_monocle3_diff_genes.pdf')
plotMarkerHeat(sc@assays$RNA@data[,order(sc$monocle3_pseudotime, decreasing = FALSE)], 
	kk$classification[order(sc$monocle3_pseudotime, decreasing = TRUE)], 
	names(sort(gene_order)), averageCells = 10^1, clusterGenes = FALSE, 
	gap = FALSE, outs = FALSE, plotheat = TRUE)
dev.off()



#--- Functional enrichment analysis of cell partitions based on pseudotimes

# Split ordered cells in 4 equal sized bins
cells = colnames(sc@assays$RNA@data[,order(sc$monocle3_pseudotime, decreasing = FALSE)])
n_bins = 4
partitions = split(cells, sort(1:length(cells) %% n_bins))
partition.length = lapply(partitions, length)
partition.numbers = rep(names(partitions), partition.length)
names(partition.numbers) = cells
sc$pseudo_partitions = partition.numbers[colnames(sc)]


# Differentially expressed genes between partitions
Idents(sc) = 'pseudo_partitions'
partition.markers = FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, 
		logfc.threshold = 0.25, test.use = 'wilcox')

data = partition.markers %>% filter(p_val_adj < 0.05) %>% 
	group_by(cluster) %>%
	arrange(desc(avg_log2FC))



# Functional enrichment analysis (only overrepresentation)
# Get top 10 terms in each group
term.list = functional_enrichment(data, organism = 'hsapiens', n = 10, colnames = c(0:3))


# Plot heatmap
pdf(file = 'fib_pseudotime_GO_heatmap.pdf')
heatmap.2(t(term.list), dendrogram = 'none', trace = 'none', scale = 'none', 
	density = 'none', Colv = FALSE, Rowv = TRUE,
	col = colorRampPalette(c('white', 'dodgerblue4'))(n = 100), key = FALSE, 
	cexRow = 0.7, cexCol = 0.7, margins = c(5,20))
dev.off()



# Pathway enrichment analysis with PID data

# Custom annotation ID is gp__cioE_CsPz_aj0 (use as organism)
gmt.path = '../data/20201118_c2.cp.pid.v7.2.symbols.gmt'
upload_GMT_file(gmtfile = gmt.path)

# Get top 10 terms in each group
term.list = functional_enrichment(data, organism = 'gp__cioE_CsPz_aj0', n = 10, colnames = c(0:3))


# Plot heatmap
pdf(file = 'fib_pseudotime_PID_heatmap.pdf')
heatmap.2(t(term.list), dendrogram = 'none', trace = 'none', scale = 'none', 
	density = 'none', Colv = FALSE, Rowv = TRUE,
	col = colorRampPalette(c('white', 'dodgerblue4'))(n = 100), key = FALSE, 
	cexRow = 0.7, cexCol = 0.7, margins = c(5,20))
dev.off()



#---- Add collagen score to UMAP embedding of fibroblasts to validate trajectory directionality
NABA = processNABA()
NABA_SETS = names(NABA)
ctrl_genes = 35
DefaultAssay(sc) = 'RNA'

for (gset in NABA_SETS){
	features = list(gset = NABA[[gset]])
	sc = AddModuleScore(object = sc, features = features, name = gset, ctrl = ctrl_genes)

	pdf(file = paste0('fib_', gset, '.pdf'))
	plot = FeaturePlot(sc, features = paste0(gset, '1'), label = FALSE, combine = FALSE)
	print(wrap_plots(plot) +
		scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red',
				midpoint = 0, limits = c(-0.5,0.5), oob = scales::squish))
	dev.off()
}



#---- Add progenitor and end-state marker genes to UMAP embedding of fibroblasts
genes = c('POSTN', 'SCARA5', 'COL1A1', 'FN1', 'COL15A1')
for (gene in genes){
	pdf(file = paste0('fib_', gene, '.pdf'))
	print(FeaturePlot(sc, feature = gene))
	dev.off()
}


