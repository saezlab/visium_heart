###{r}
library(Matrix)
library(Seurat)
library(optparse)
library(stringr)
library(dplyr)
library(RANN)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=100)

suppressPackageStartupMessages(library(future))
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 100000 * 1024^2)


library(Rcpp)


Rcpp::sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);    } 

    return(cor);

  }'
)

AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-s", "--SAMPLE"), type="character", default="CK166",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)

atac_name  = pa$SAMPLE

#rna_name = paste0("CK", as.numberic(substr(atac_name, 3, 5))-8)
atac_to_rna = setNames( c("CK158", "CK159", "CK160",  "CK162", "CK163", "CK164","CK165"), 
                    c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174"))

rna_name = atac_to_rna[atac_name] 


message("Scaling rna", date())

rna_names <- sprintf("save/aggred_counterpart_%s.Rds", rna_name)
rna_mtx <- readRDS(file=rna_names)

message("rna ", date())
rna_seurat <- CreateSeuratObject(counts = rna_mtx)
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)
rna_seurat <- ScaleData(rna_seurat, features=rownames(rna_seurat))
scaled_rna_mtx <- as.matrix(rna_seurat@assays$RNA@scale.data)


message("Scaling atac", date())
atac_names <- sprintf("save/%s_aggred_atac.Rds", atac_name)
atac_mtx <- readRDS(file=atac_names)
atac_seurat <- CreateSeuratObject(counts = atac_mtx)
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- FindVariableFeatures(atac_seurat)
atac_seurat <- ScaleData(atac_seurat, features=rownames(atac_seurat))
scaled_atac_mtx <- as.matrix(atac_seurat@assays$RNA@scale.data)


saveRDS(scaled_rna_mtx, file=sprintf("save/%s_scale_rna_mtx.Rds", rna_name))
saveRDS(scaled_atac_mtx, file=sprintf("save/%s_scale_atac_mtx.Rds", atac_name))

#scaled_rna_mtx <- readRDS(file=sprintf("save/%s_scale_rna_mtx.Rds", rna_name))
#scaled_atac_mtx <- readRDS(file=sprintf("save/%s_scale_atac_mtx.Rds", atac_name))



new_ann <- data.frame(peak=character(), gene=character(), 
                      distance=numeric(), #peak_type=character(), 
                      corr=numeric(), stringsAsFactors = F)


ann <- read.csv(file=sprintf("save/peak2gene_putative-simple_%s.tsv", atac_name), sep="\t", stringsAsFactors = F)
###

###{r}
#future_sapply()
#scaled_atac_mtx <- as.matrix(atac_mtx)
#scaled_rna_mtx <- as.matrix(rna_mtx)

#rna_genes <-rownames(scaled_rna_mtx)
#names(rna_genes) <- rna_genes


`%notin%` <- Negate(`%in%`) 

atac_rownames <- rownames(scaled_atac_mtx)
atac_rownames <- sub("-", ":", atac_rownames)
rownames(scaled_atac_mtx) <- atac_rownames

ann <- ann %>% filter(gene %in% rownames(scaled_rna_mtx))
ann <- ann %>% filter(peak %in% rownames(scaled_atac_mtx))


peaks <- ann$peak
genes <- ann$gene

peak_idx <- match(peaks, rownames(scaled_atac_mtx))
gene_idx <- match(genes, rownames(scaled_rna_mtx))


corrs <- rowCorCpp(peak_idx, gene_idx, scaled_atac_mtx, scaled_rna_mtx)
new_ann <- data.frame(peaks=ann$peak, gene=ann$gene, distance=ann$distance, corr=corrs, stringsAsFactors=F)


order_v <- factor(sapply(new_ann[["peaks"]], function(x) str_split(x, ":")[[1]][1] ), 
                    levels=c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
                             'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 
                             'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                             'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'))


new_ann <- new_ann[order(order_v) ,]



message( paste0(atac_name,  "  finished calculate: " , dim(ann)[1], "  ", date()))
fname = paste0("save/corr_", atac_name, ".tsv")
write.table(new_ann, file = fname, sep="\t",  quote=FALSE)
###
