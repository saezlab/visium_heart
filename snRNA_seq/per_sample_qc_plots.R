# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Generate per sample QC plots

library(Seurat)
library(rlist)
library(cowplot)
library(ggplot2)
samples = c('CK158', 'CK159', 'CK160', 'CK162', 'CK163', 'CK164', 'CK165')



obj.list = c()
for (sample in samples){
	sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
	obj.list = list.append(obj.list, sc)
}
sc.obj = merge(obj.list[[1]], obj.list[-1], add.cell.ids = samples)



# Plot cell count per sample in bar chart
cell.count = as.data.frame(table(sc.obj$orig.ident))
colnames(cell.count) = c('Sample', 'Count')

p1 = ggplot(data = cell.count, aes(x = Sample, y = log10(Count))) +
    geom_bar(aes(color = Sample, fill = Sample), 
             alpha = 1, stat = 'identity') +
    xlab('') + ylab('Number of \nValid Cells (log10)') +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none')


# Plot UMI and gene count per sample in violin plot
umi = data.frame(cell = names(Idents(sc.obj)), 
	UMI = as.numeric(sc.obj$nCount_RNA), gene = as.numeric(sc.obj$nFeature_RNA),
	sample = sc.obj$orig.ident)

p2 = ggplot(data = umi, aes(x = sample, y = log10(UMI))) +
    geom_boxplot(aes(color = sample)) +
    geom_violin(aes(color = sample, fill = sample), alpha = 0.5) +
    xlab('') + ylab('Number of \nUMIs (log10)') +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none')


p3 = ggplot(data = umi, aes(x = sample, y = log10(gene))) +
    geom_boxplot(aes(color = sample)) +
    geom_violin(aes(color = sample, fill = sample), alpha = 0.5) +
    xlab('') + ylab('Number of \nGenes (log10)') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          axis.ticks.x = element_blank(),
          legend.position = 'none')



pdf(file = 'scRNA_QC.pdf', width = 4, height = 6)
plot_grid(p1, p2, p3, ncol = 1, align = 'v', rel_heights = c(1, 1, 1.3))
dev.off()


