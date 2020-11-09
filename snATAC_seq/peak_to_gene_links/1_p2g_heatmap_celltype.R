#rm(list=ls())
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Rcpp))
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
#plan("multiprocess", workers = 20) 
#options(future.globals.maxSize = 100000 * 1024^2)


##--------------------------------------------------------------------------------------
## changed from https://github.com/AEBilgrau/correlateR/blob/master/src/corFamily.cpp
##--------------------------------------------------------------------------------------

Rcpp::sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // [[Rcpp::export]]
  Rcpp::NumericMatrix corRcpp(Rcpp::NumericMatrix & X) {
  
  const int m = X.ncol();
  const int n = X.nrow();
  
  // Centering the matrix
  //X = centerNumericMatrix(X);
  for(int j=0; j < m; ++j){
      X(Rcpp::_, j ) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }

  
  Rcpp::NumericMatrix cor(m, m);
  
  // Degenerate case
  if (n == 0) {
    std::fill(cor.begin(), cor.end(), Rcpp::NumericVector::get_na());
    return cor; 
  }
  
  // Compute 1 over the sample standard deviation
  Rcpp::NumericVector inv_sqrt_ss(m);
  for (int i = 0; i < m; ++i) {
    inv_sqrt_ss(i) = 1/sqrt(Rcpp::sum(X(Rcpp::_, i)*X(Rcpp::_, i)));
  }
  
  // Computing the correlation matrix
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      cor(i, j) = Rcpp::sum(X(Rcpp::_,i)*X(Rcpp::_,j)) *
        inv_sqrt_ss(i) * inv_sqrt_ss(j);
      cor(j, i) = cor(i, j);
    }
  }
  return cor;
}')





parser <- OptionParser()

parser <- add_option(parser, c("-s", "--sample"), type="character", default="CK166",
                     help="patient code(HD, Health, HTN)  [default %default]",
                     metavar="CHAR")

parser <- add_option(parser, c("-n", "--number"), type="character", default="9",
                     help="Row cluster number [default %default]",
                     metavar="INTEGER")

parser <- add_option(parser, c("-g", "--generate"), type="character", default="NO",
                    help="Generate seurat clusters [default %default]",
                    metavar="character")

parser <- add_option(parser, c("-d", "--dist"), type="character", default="NO",
                    help="calculate df.dist from correlation [default %default]",
                    metavar="character")

parser <- add_option(parser, c("-p", "--plots"), type="character", default="",
                    help="plots for different pval_option [default %default]",
                    metavar="character")

parser <- add_option(parser, c("-r", "--rowclusters"), type="character", default="NO",
                    help="use the calculated rowclusters or not [default %default]",
                    metavar="character")

pa = parse_args(parser)

number.clusters <- as.numeric(pa$number)

seurat.generate <- pa$generate
#read_from_rds <- pa$RDS
atac_name = pa$type 
calc_dist <- pa$dist
row_clusters <- pa$rowclusters
atac_name <- pa$sample

num_dict <- c(
           "CK166" = 8,
           "CK167" = 8,
           "CK168" = 8,
           "CK170" = 6,
           "CK171" = 9,
           "CK173" = 10,
           "CK174" = 6 
    )

number.clusters <- num_dict[atac_name]



plots_opt <-  ifelse(pa$plots == "",  "",  paste0("_", pa$plots))
print(plots_opt)
dir.create(sprintf("plots%s", plots_opt))
print(sprintf("plots%s", plots_opt))


##----------------- this is for temporary --------------------------------
#plots_opt  <- "_0.01_0.5"


atac_mtx <- readRDS(file=sprintf("save/%s_scale_atac_mtx.Rds", atac_name))

atac_rownames <- rownames(atac_mtx)
atac_rownames <- sub("-", ":", atac_rownames)
rownames(atac_mtx) <- atac_rownames
#rna_mtx <- readRDS(file=sprintf("%s_scale_rna_mtx", atac_name))


fdf <- read.csv(file=sprintf("save/corr_pval_final_%s%s.tsv",atac_name, plots_opt),
                         sep="\t", stringsAsFactors = F)

order_atac_mtx <- atac_mtx[fdf$peak, ]#[1:1000, 1:1000]
rm(atac_mtx)
gc()
#
####-------------------celltype order begin---------------------------
atac_dir <- paste0("../ATAC_Single_Sample_V2/data/", atac_name)
atac <- readRDS(file=file.path(atac_dir, paste0(atac_name, "_unionPeaks.Rds")))



meta_data <- atac@meta.data
meta_data <- meta_data[which(rownames(meta_data) %in% colnames(order_atac_mtx)), ]
meta_data <- meta_data[order(meta_data$celltype), ]
lst <-  list()

#tbl <- table(meta_data$celltype)
#nms <- names(tbl)
meta_data$celltype <- as.character(meta_data$celltype)
ocelltypes <- sort(unique(meta_data$celltype))


for (i in 1:length(ocelltypes)){
    celltype <- ocelltypes[i]
    cells <- rownames(meta_data[meta_data$celltype == celltype, ])
    lst[[i]] <- cells    
}

ordered_lst = list()
for(i in 1:length(lst)){
    a_mtx <- order_atac_mtx[, lst[[i]]]
    corr <- corRcpp(a_mtx)
    c_dist <- as.dist(corr)
    c_clusters <- hclust(1 - c_dist)
    c_order <- c_clusters$order 
    ordered_lst[[i]] <- lst[[i]][c_order]
}

celltype_column_order <- unlist(ordered_lst)
###----------------annotation begin--------------------------

cn <- celltype_column_order
celltypes <- meta_data[cn,]$celltype

ha = HeatmapAnnotation(celltype = celltypes) #,
#    col = list(celltype = c("adipocytes" = "#8dd3c7", "cardiomyocytes" = "#bebada", "endothelium" = "#fb8072","fibroblasts"="#80b1d3", "lymphatic_endothelium"="#fdb462", "mast_cells"="#b3de69", "monocytes"="#fccde5", "neuronal"="#d9d9d9", "pericytes"="#bc80bd", "T-cells"="#ccebc5", "VSMCs"="#ffed6f")))

###----------------annotation ends--------------------------


colors <- c("#f7fcf0",  "#e0f3db",  "#ccebc5",  "#a8ddb5",  "#7bccc4",  "#4eb3d3",  "#2b8cbe",  "#0868ac",  "#084081")
col_fun = colorRamp2(c(-4:4), colors)

#order_atac_mtx <- as.matrix(sobj@assays$RNA@scale.data)
#rownames(order_atac_mtx) #<- gsub("-", "_", rownames(order_atac_mtx))

if(calc_dist == "YES"){
    print(paste(atac_name, ": dist calculating", date()))
    rns <- rownames(order_atac_mtx)
    corr_c_version = corRcpp(t(order_atac_mtx)) 
    rownames(corr_c_version) <- rns 
    colnames(corr_c_version) <- rns 
    df.dist <- 1 - corr_c_version

    saveRDS(df.dist, file=sprintf("plots%s/%s_df.dist.Rds", plots_opt, atac_name))
    print(paste(atac_name, ": dist calculated", date()))
}else{
    print(paste(atac_name, ":loading dist", date()))
    df.dist = readRDS(file=sprintf("plots%s/%s_df.dist.Rds", plots_opt, atac_name))
}


if(row_clusters == "YES"){
    print(paste(atac_name, ": pamming", date()))
    pa = pam(df.dist, k = number.clusters, diss = TRUE, cluster.only=TRUE)
    rm(df.dist)
    gc()
    print(paste(atac_name, ": pammed", date()))
    row_split = pa
    saveRDS(row_split, file=sprintf("plots%s/%s_pa_%d.Rds", plots_opt,atac_name, number.clusters))
}else{
    print(paste(atac_name, ":loading rowclusters", date()))
    row_split <- readRDS(file=sprintf("plots%s/%s_pa_%d.Rds", plots_opt,atac_name, number.clusters)) 
}

pdf(sprintf("plots%s/%s_heatmap_celltype_%d.pdf", 
            plots_opt, atac_name, number.clusters), width = 8, height = 8)

print(paste(atac_name, ": heatmapping", date()))
h1 <- Heatmap(order_atac_mtx[, celltype_column_order], 
          split = row_split,
          cluster_rows = T, show_row_dend = F,
          #column_split = sobj$seurat_clusters,
          show_column_dend = F, cluster_columns = F,
          show_row_names = F,  
          show_column_names = F,
          name = sprintf("ATAC %s", atac_name),
          col = col_fun,
          use_raster=T,
          raster_device="png",
          top_annotation = ha
        )

draw(h1)

dev.off()
print(paste(atac_name, ": heatmapped", date()))

#col_order <- column_order(h1)
row_order <- row_order(h1)


saveRDS(ordered_lst, file=sprintf("plots%s/col_order_celltype_%d_%s.Rds", 
                     plots_opt, number.clusters, atac_name))
saveRDS(row_order, file=sprintf("plots%s/row_order_celltype_%d_%s.Rds",
                     plots_opt, number.clusters, atac_name))

