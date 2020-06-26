suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr)) ## for str_remove
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(fclust))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(tiff))

plan("multiprocess", workers = 20)
options(future.globals.maxSize = 80000 * 1024^2)
#


#----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

parser <- OptionParser()

parser <- add_option(parser, c("-s", "--sample"), type="character", default="CK174",
                     help="patient atac name [default %default]",
                     metavar="CHAR")
parser <- add_option(parser, c("-c", "--color"), type="character", default="A",
                     help="colors option [default %default]",
                     metavar="CHAR")
parser <- add_option(parser, c("-p", "--plots"), type="character", default="",
                    help="plots for different pval_option [default %default]",
                    metavar="character")


cols = c('Cardiomyocytes' = '#800000',
'Cardiomyocytes 1' = '#800000',
'Cardiomyocytes 2' = '#9A6324',
'Cardiomyocytes 3' = '#808000',
'Fibroblasts' = '#911eb4',
'Fibroblasts 1 COL15A1+' = '#911eb4',
'Fibroblasts 1' = '#911eb4',
'Fibroblasts 2 SCARA5+' = '#e6beff',
'Fibroblasts 2' = '#e6beff',
'Fibroblasts 3' = '#f032e6',
'Endothelial cells' = '#000075',
'Endothelial cells 1' = '#000075',
'Endothelial cells 2' = 'blue',
'Endothelial cells 2 POSTN+' = 'blue',
'Endothelial cells 3' = '#568198',
'Endothelial cells 3 PLVAP+' = '#568198',
'Endothelial cells 3 VEGFC+' = '#568198',
'Endothelial cells 4' = '#469990',
'Endothelial cells 4 SEMA3G+' = '#469990',
'Macrophages' = '#e6194B',
'Macrophages 1 CD163+' = '#e6194B',
'Macrophages 2 CD11C+' = '#fabebe',
'Pericytes' = '#f58231',
'Pericytes EGFLAM+' = '#f58231',
'T cells' = '#ffe119',
'Lymphatic endothelial cells' = '#ffd8b1',
'Adipocytes' = '#000000',
'Neuronal cells' = '#42d4f4',
'Erythrocytes' = '#999999',
'Proliferating cells' = '#999999',
'Damaged endothelial cells' = '#999999',
'Vascular smooth muscle cells' = '#aaffc3')


pa = parse_args(parser)

color.option <- pa$color 
atac_name = pa$sample

atac_to_rna = setNames( c("CK158", "CK159", "CK160",  "CK162", "CK163", "CK164","CK165"), 
                    c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174"))
rna_name = atac_to_rna[atac_name]

plots_opt <-  ifelse(pa$plots == "",  "",  paste0("_", pa$plots))


print(plots_opt)
dir.create(sprintf("plots%s", plots_opt))
print(sprintf("plots%s", plots_opt))



atac_mtx <- readRDS(file=sprintf("save/%s_scale_atac_mtx.Rds", atac_name))

atac_rownames <- rownames(atac_mtx)
atac_rownames <- sub("-", ":", atac_rownames)
rownames(atac_mtx) <- atac_rownames

rna_mtx <- readRDS(file=sprintf("save/%s_scale_rna_mtx.Rds", rna_name))

fdf <- read.csv(file=sprintf("save/corr_pval_final_%s%s.tsv", atac_name, plots_opt),
                                         sep="\t", stringsAsFactors = F)

order_atac_mtx <- as.matrix(atac_mtx[fdf$peak, ])#[1:1000, 1:1000]



if(plots_opt == ""){
cluster_orders <- list(
"CK166" = c(8,1,2,6,3,7,5,4),
"CK167" = c(5,3,8,2,1,4,7,6),
"CK168" = c(8,6,5,1,2,3,4,7),
"CK170" = c(5,2,4,1,3,6),
"CK171" = c(1,8,6,7,5,4,9,2,3),
"CK173" = c(3,2,1,6,5,8,9,7,10,4),
"CK174" = c(2,6,3,4,5,1)
) ### pval=0.05 
}
if(plots_opt == "_0.01"){
cluster_orders <- list(
    "CK166" = c(1,5,4,8,2,7,3,6),
    "CK167" = c(4,8,7,2,1,3,6,5),
    "CK168" = c(8,5,3,1,2,4,6,7),
    "CK170" = c(3,1,2,5,4,6),
    "CK171" = c(1,7,4,3,5,6,2,8,9),
    "CK173" = c(1,2,9,8,4,7,6,5,10,3),
    "CK174" = c(6,2,3,5,4,1)
    )   ##pval0.01
}

if(plots_opt == "_0.01_0.3"){
cluster_orders <- list(
    "CK166" = c(4,3,8,7,1,2,5,6),
    "CK167" = c(5,1,7,3,2,6,8,4),
    "CK168" = c(8,5,3,1,2,6,4,7),
    "CK170" = c(1,3,4,6,2,5),
    "CK171" = c(2,5,9,1,6,7,4,3,8),
    "CK173" = c(1,7,8,2,4,6,5,10,9,3),
    "CK174" = c(2,4,3,6,5,1) 
    )  ## corr0.3

}
if(plots_opt == "_0.01_0.5"){
cluster_orders <- list(
    "CK166" = c(8,4,2,6,1,3,7,5,9))   ##0.5 
}


number.clusters = length(cluster_orders[[atac_name]])

row_order_list <- readRDS(sprintf("plots%s/row_order_celltype_%d_%s.Rds", plots_opt,number.clusters, atac_name))

row_order_ <- unlist(row_order_list[as.character(1:9)]) 



rorder <- as.character(cluster_orders[[atac_name]])
nlist <- row_order_list[rorder]
row_order_ <- unlist(nlist)
row_order_anno <-lapply(seq_along(nlist),
                function(y, n, i) { rep(n[[i]] ,length(y[[i]])) },
                y=nlist, n=names(nlist))

number.clusters <- length(row_order_list)



col_order_list <- readRDS(sprintf("plots%s/col_order_celltype_%d_%s.Rds", plots_opt, number.clusters, atac_name))

col_order_ <- unlist(col_order_list)


#colors <- c("#f7fcf0",  "#e0f3db",  "#ccebc5",  "#a8ddb5",  "#7bccc4",  "#4eb3d3",  "#2b8cbe",  "#0868ac",  "#084081")
#col_fun = colorRamp2(c(-4:4), colors)


cn <-unlist(col_order_list)

atac_dir <- paste0("../ATAC_Single_Sample_V2/data/", atac_name)
atac <- readRDS(file=file.path(atac_dir, paste0(atac_name, "_unionPeaks.Rds")))
meta_data <- atac@meta.data
celltypes <- as.character(meta_data[cn,]$celltype)

table(celltypes)

ha = HeatmapAnnotation(celltype = celltypes,
                         annotation_name_gp = gpar(fontsize = 0),
                         col = list(celltype = cols))

lens <- sapply(row_order_anno, length)
x = 1
vec = c(x)
for (i in 1:(length(lens)-1)){
    x = x + lens[i]
    vec = c(vec, x)
}

if(plots_opt == ""){

    cluster_anno <- switch(atac_name,
"CK166" = c("8" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "2" = "Endothelial cells 1",
  "6" = "Endothelial cells 2 POSTN+",
  "3" = "Fibroblasts 1 COL15A1+",
  "7" = "Macrophages 1 CD163+",
  "5" = "Pericytes",
  "4" = "Pericytes"),
  
 "CK167" = c("5" = "Cardiomyocytes",
  "3" = "Endothelial cells 1",
  "8" = "Endothelial cells 1",
  "2" = "Endothelial cells 2 POSTN+",
  "1" = "Fibroblasts",
  "4" = "Fibroblasts",
  "7" = "Macrophages",
  "6" = "Pericytes"),


"CK168" = c("8" = "Cardiomyocytes 1",
  "6" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "3" = "Macrophages",
  "4" = "Macrophages",
  "7" = "Pericytes EGFLAM+"),

"CK170" = c("5" = "Cardiomyocytes",
  "2"= "Cardiomyocytes",
  "4" = "Damaged endothelial cells",
  "1" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "6" = "Fibroblasts"),

"CK171" = c("1" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "6" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "5" = "Fibroblasts",
  "4" = "Macrophages",
  "9" = "Neuronal cells",
  "2" = "Pericytes",
  "3" = "Vascular smooth muscle cells"),

"CK173" = c("3" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "6" = "Endothelial cells 1",
  "5" = "Endothelial cells 1",
  "8" = "Fibroblasts 1",
  "9" = "Fibroblasts 1",
  "7" = "Fibroblasts 2 SCARA5+",
  "10" = "Macrophages",
  "4" = "Pericytes"),

"CK174" = c("2" = "Cardiomyocytes 1",
  "6" = "Cardiomyocytes 2",
  "3" = "Endothelial cells",
  "4" = "Fibroblasts",
  "5" = "Macrophages",
  "1" = "Pericytes"))
}

if(plots_opt == "_0.01"){
    # 4,6,5,9,3,1,2,8,7

   cluster_anno <- switch(atac_name,
   "CK166" = c("1" = "Cardiomyocytes",
  "5" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "8" = "Endothelial cells 2 POSTN+",
  "2" = "Fibroblasts 1 COL15A1+",
  "7" = "Macrophages 1 CD163+",
  "3" = "Pericytes",
  "6" = "Pericytes"),
  
 "CK167" = c("4" = "Cardiomyocytes",
  "8" = "Endothelial cells 1",
  "7" = "Endothelial cells 2 POSTN+",
  "2" = "Endothelial cells 3 PLVAP+",
  "1" = "Fibroblasts",
  "3" = "Fibroblasts",
  "6" = "Macrophages",
  "5" = "Pericytes"),

"CK168" = c("8" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 1",
  "3" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "4" = "Fibroblasts 1",
  "6" = "Macrophages",
  "7" = "Pericytes EGFLAM+"),

"CK170" = c("3" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "2" = "Damaged endothelial cells",
  "5" = "Endothelial cells 1",
  "4" = "Fibroblasts",
  "6" = "Fibroblasts"),

"CK171" = c("1" = "Cardiomyocytes",
  "7" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "5" = "Macrophages",
  "6" = "Neuronal cells",
  "2" = "Pericytes",
  "8" = "Pericytes",
  "9" = "Vascular smooth muscle cells"),

"CK173" = c("1" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "9" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "6" = "Fibroblasts 1",
  "5" = "Fibroblasts 2 SCARA5+",
  "10" = "Macrophages",
  "3" = "Pericytes"),

"CK174" = c("6" = "Cardiomyocytes 1",
  "2" = "Endothelial cells",
  "3" = "Fibroblasts",
  "5" = "Fibroblasts",
  "4" = "Macrophages",
  "1" = "Pericytes")
        )

    
}

if(plots_opt == "_0.01_0.3"){
    #3,1,4,6,8,2,7,9,5

cluster_anno <- switch(atac_name,
"CK166" = c("4" = "Cardiomyocytes",
  "3" = "Endothelial cells 1",
  "8" = "Endothelial cells 1",
  "7" = "Endothelial cells 2 POSTN+",
  "1" = "Fibroblasts 1 COL15A1+",
  "2" = "Fibroblasts 1 COL15A1+",
  "5" = "Macrophages 1 CD163+",
  "6" = "Pericytes"),
  
  
  
"CK167" = c("5" = "Cardiomyocytes",
  "1" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "2" = "Fibroblasts",
  "6" = "Macrophages",
  "8" = "Neuronal cells",
  "4" = "Pericytes"),


"CK168" = c("8" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 1",
  "3" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "6" = "Macrophages",
  "4" = "Pericytes EGFLAM+",
  "7" = "Pericytes EGFLAM+"),


"CK170" = c("1" = "Cardiomyocytes",
  "3" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "6" = "Endothelial cells 1",
  "2" = "Fibroblasts",
  "5" = "Fibroblasts"),



"CK171" = c("2" = "Cardiomyocytes",
  "5" = "Endothelial cells 1",
  "9" = "Endothelial cells 1",
  "1" = "Fibroblasts",
  "6" = "Macrophages",
  "7" = "Neuronal cells",
  "4" = "Pericytes",
  "3" = "Pericytes",
  "8" = "Vascular smooth muscle cells"),



"CK173" = c("1" = "Cardiomyocytes",
  "7" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "6" = "Endothelial cells 1",
  "5" = "Fibroblasts 1",
  "10" = "Fibroblasts 1",
  "9" = "Macrophages",
  "3" = "Pericytes"),


"CK174" = c("2" = "Cardiomyocytes 1",
  "4" = "Cardiomyocytes 2",
  "3" = "Endothelial cells",
  "6" = "Fibroblasts",
  "5" = "Macrophages",
  "1" = "Pericytes")

)
}
if(plots_opt == "_0.01_0.5"){
    cluster_anno <- c(
    "8" = "Cardiomyocytes",
    "4" = "Fibroblasts 1 COL15A1+",
    '2' = "Macrophages 1 CD163+",
    '6' = "Endothelial cells 1",
    '1' = "Proliferating cells",
    '3' = "Pericytes",
    '7' = "Pericytes",
    '5' = "Endothelial cells 2 POSTN+",
    '9' = "Neuronal cells"
)
}
#5,9,6,7,4,2,3,8,1





annotated_cluster <- cluster_anno[unlist(row_order_anno)]

ca = rowAnnotation(clusters = annotated_cluster,
                   annotation_name_gp = gpar(fontsize = 0),
                   col = list(clusters = cols),
                   show_legend = FALSE
                   )



all_peaks <- rownames(order_atac_mtx)[row_order_]


col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))

pdf(sprintf("plots%s/reorder_%s_%s_heatmap_%d.pdf", plots_opt ,color.option,atac_name , number.clusters), width = 8, height = 8)

h1 <- Heatmap(order_atac_mtx[all_peaks, col_order_],
          cluster_rows = F, show_row_dend = F,
          #split = cluster_split,
          #row_order=r_o,
          #row_title = str_c(rev(rorder), collapse=" <- "),
          show_column_dend = F, cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          column_title = "ATAC",
          name = sprintf("ATAC %s", atac_name),
          top_annotation = ha,
          col = col_fun,
          left_annotation = ca,
          use_raster=T,
          raster_device="CairoTIFF"
        )



#draw(h1)

rna_dir <- "../scRNA_filtered"
rna <- readRDS(file=file.path(rna_dir, paste0(rna_name, ".filtered.annotated.rds")))

rna$celltype <- Idents(rna)

col_order_i <- match(col_order_, colnames(order_atac_mtx))
rna_metadata <- rna@meta.data

cn_rna <- colnames(rna_mtx)
cn_rna <- str_remove(cn_rna, "\\.\\d+$")
rna_celltypes  <- rna_metadata[cn_rna,]$celltype
col_fun <-  circlize::colorRamp2(c(-2, 0, +2), c("purple", "black", "yellow"))
ha = HeatmapAnnotation(rna_celltype = rna_celltypes[col_order_i],
                         show_legend = F,
                         annotation_name_gp = gpar(fontsize = 0),
                         col = list(rna_celltype = cols))

h2 <- Heatmap(rna_mtx[fdf$gene[row_order_], col_order_i],
          cluster_rows = F, show_row_dend = F,
          cluster_columns = F, show_column_dend = F,
          show_row_names = F,
          column_title = "RNA",
          show_column_names = F,
          name = sprintf("RNA %s", atac_name),
          top_annotation=ha,
          col = col_fun,
          use_raster=T,
          raster_device="CairoTIFF"
        )

draw(h1 + h2)
dev.off()


peaks_clusters <- lapply(row_order_list, function(x)rownames(order_atac_mtx)[x] )
genes_clusters <- lapply(row_order_list, function(x) rownames(rna_mtx[fdf$gene, ])[x]) 

merge_info <- Map(list, peaks_clusters, genes_clusters)
fname_merge_info <- sprintf("plots%s/merge_info_%s_%d.Rds", plots_opt,atac_name, number.clusters)
saveRDS(merge_info, file=fname_merge_info)



fname_atac <- sprintf("plots%s/final_mtx_atac_%s_%d.Rds", plots_opt,atac_name, number.clusters)
fname_rna <- sprintf("plots%s/final_mtx_rna_%s_%d.Rds", plots_opt, rna_name, number.clusters)
saveRDS(order_atac_mtx[all_peaks, col_order_], file=fname_atac)
saveRDS(rna_mtx[fdf$gene[row_order_], col_order_i], file=fname_rna)

