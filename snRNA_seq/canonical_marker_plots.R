# Copyright (c) [2020] [Monica T. Hannani]
# mhannani@ukaachen.de


#---- Violin plots of canonical cell type markers per sample
library(Seurat)
library(ggplot2)
library(cowplot)
source('./sc_source/sc_source.R')



get_violin = function(object, features.use){
    p = VlnPlot(object = object, 
                 features = features.use, 
                 pt.size = 0, cols = cell.type.colors,
                 combine = FALSE)
    
    ps = lapply(p, function(x) x + coord_flip() + NoLegend() +
                     theme_bw() +
                     theme(plot.title = element_text(angle = 90), legend.position = 'none',
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           plot.margin = unit(c(0, 0, 0, 0), 'cm'),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_rect(colour = 'black', size = 1)))
    
    p = plot_grid(plotlist = ps, nrow = 1, align = 'h')
    
    return(p)
}



# CK158 (P1 C)
sample = 'CK158'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('MEG3', 'COL15A1', 'NEGR1',
                'MRC1', 'SIGLEC1',
                'VWF', 'TPO', 'FLT1',
                'MYOM2', 'RBM20', 'TRDN',
                'SCARA5', 'PCOLCE2',
                'EGFLAM', 'ABCC9', 'CARMN',
                'POSTN', 'NRG3', 'EMCN',
                'NOTCH4', 'VEGFC',
                'ITGAX', 'ZNF804A', 'RUNX1',
                'SKAP1', 'CD247', 'ITK',
                'NRXN1', 'NRXN3', 'GRIK3',
                'MYH11', 'SLC38A11', 'LMOD1',
                'ASPM', 'MKI67')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK159 (P5 CZ)
sample = 'CK159'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('COL15A1', 'COL3A1', 'MEG3',
                'VWF', 'EMCN', 'BTNL9',
                'EGFLAM', 'RGS5', 'PDGFRB', 'ABCC9',
                'MRC1', 'CD163', 'MERTK',
                'RBM20', 'LDB3', 'MYO18B',
                'POSTN', 'NOSTRIN',
                'PLVAP',
                'SEMA3G', 'NOTCH4', 'VEGFC',
                'NRXN1', 'NEGR1',
                'PTPRC', 'SKAP1', 'BCL11B')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK160 (P3 BZ)
sample = 'CK160'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('COL3A1', 'RUNX1', 'COL1A2','TNNT2',
                'PTPRC', 'CD163',
                'VWF', 'FLT1', 'EMCN', 'EGFL7',
                'PCOLCE2', 'SCARA5',
                'MYOM1', 'LDB3',
                'NOTCH3', 'RGS5', 'PDGFRB',
                'CD247', 'LCK', 'THEMIS',
                'RAMP3', 'TPO',
                'NRXN1', 'ZNF536',
                'MYH11', 'ITGA8',
                'CASC5', 'TPX2',
                'MMRN1', 'TBX1', 'FLT4', 'PROX1')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK162 (P4 CZ)
sample = 'CK162'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('CD36', 'EMCN', 'VWF', 'TPO',
                'DCN', 'COL3A1',
                'HSPH1', 'HSPD1',
                'VEGFC', 'PCSK5',
                'FOXP1', 'RAMP3',
                'ANGPT2', 'CMIP',
                'RGS5', 'EGFLAM', 'ABCC9', 'NOTCH3',
                'MYH11', 'LMOD1',
                'NEGR1', 'DCN',
                'RBM20', 'RYR1', 'MYOM1',
                'CLMP', 'GLIS3', 'RUNX1',
                'CD163',
                'PTPRC', 'SKAP1', 'ITK',
                'ZNF536', 'NRXN1',
                'MMRN1', 'PROX1')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK163 (P2 BZ)
sample = 'CK163'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('NEGR1', 'DCN', 'GLIS3',
                'VWF', 'EMCN', 'FLT1',
                'EGFLAM', 'RGS5', 'PDGFRB', 'ABCC9',
                'RBM20', 'MYO18B', 'MYOM2', 'TRDN',
                'CD163', 'MRC1', 'SIGLEC1',
                'EGFL7',
                'RAMP3', 'TPO', 'PTPRB',
                'MYH11', 'CARMN', 'TBX2',
                'SKAP1', 'CD247',
                'NRXN1', 'NRXN3', 'ZNF536',
                'SLC4A1', 'FHDC1', 'HBA1', 'HBA2')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK164 (P3 RZ)
sample = 'CK164'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('GLIS3', 'SLIT2',
                'TNNT2', 'MYO18B',
                'CD163',
                'PCOLCE2', 'SCARA5',
                'VWF', 'FLT1', 'EMCN',
                'EGFLAM', 'NOTCH3', 'RGS5',
                'SKAP1', 'CD247',
                'POSTN', 'SLCO2A1', 'ITGA8',
                'MMRN1', 'PROX1', 'FLT4',
                'NRXN1', 'NRXN3', 'ZNF536',
                'PLIN4', 'ADIPOQ')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')



# CK165 (P2 IZ)
sample = 'CK165'
sc = readRDS(file = paste0('../data/', sample, '.filtered.annotated.rds'))
Idents(sc) = factor(Idents(sc), levels = cell.type.order)

marker.genes = c('TNNT2',
                'MYO18B',
                'VWF', 'PECAM1',
                'MEG3', 'GLI2', 'MEG8',
                'MYH11', 'NOTCH3', 'PDGFRB',
                'RBM47', 'HMOX1', 'RXRA', 'CD163')

p = get_violin(object = sc, features.use = marker.genes)
ggsave(filename = paste0(sample, '.pdf'), plot = p, width = 12, height = 6, units = 'in')

