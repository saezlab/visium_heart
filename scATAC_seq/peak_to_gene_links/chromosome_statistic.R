library(Matrix)
library(Seurat)
library(optparse)
library(RANN)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel(cores=24)


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

atac_name   = pa$SAMPLE

#rna_name = paste0("CK", as.numberic(substr(atac_name, 3, 5))-8)
atac_to_rna = setNames( c("CK158", "CK159", "CK160",  "CK162", "CK163", "CK164","CK165"), 
                    c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174"))

rna_name = atac_to_rna[atac_name]


atac_mtx_undf <- readRDS(file=sprintf("save/%s_scale_atac_mtx.Rds", atac_name))
rna_mtx_undf <- readRDS(file=sprintf("save/%s_scale_rna_mtx.Rds", rna_name))

atac_rownames <- rownames(atac_mtx_undf)
atac_rownames <- sub("-", ":", atac_rownames)
rownames(atac_mtx_undf) <- atac_rownames


#atac_mtx <- readRDS(file=sprintf("save/%s_aggred_atac.Rds", type_code))
#rna_mtx <- readRDS(file=sprintf("save/%s_reorder_rna_mtx.Rds", type_code))
#


corr_df <- read.csv(sprintf("save/corr_%s.tsv", atac_name), sep="\t", stringsAsFactors=F)
chromosome <- read.csv(sprintf("save/chromosome_range_%s.txt", atac_name), 
                                    sep = " ", stringsAsFactors = F)
nr <- nrow(chromosome)
maxidx <- chromosome$end[nr]
whole_seq = seq(1, maxidx)


chr_all_corrs <- foreach(i = 1:nrow(chromosome)) %dopar% {
  set.seed(1)
  a_chr <- chromosome[i, ]
  tseq <- setdiff(whole_seq, a_chr$start:a_chr$end)
  peaks_idx <- sample(tseq)[1:1000]
  peaks <- corr_df[peaks_idx, ]$peak
  
  genes <- strsplit(a_chr$genes, ",")[[1]]
  genes <- intersect(genes, rownames(rna_mtx_undf))
    

  peak_mtx_idx = match(peaks, rownames(atac_mtx_undf)) 
  gene_mtx_idx = match(genes, rownames(rna_mtx_undf))
  
  grid <- expand.grid(peak_mtx_idx, gene_mtx_idx)
  all_corrs <- rowCorCpp(as.integer(grid[,1]), as.integer(grid[,2]), atac_mtx_undf, rna_mtx_undf)

  all_corrs = as.vector(all_corrs)
  all_corrs = na.omit(all_corrs)
}


names(chr_all_corrs) <-  chromosome$chromosome

message( paste0(atac_name,  "  finished calculate: ", date()))
fname = paste0("save/", atac_name, "_null_hypothesis.Rds")
saveRDS(chr_all_corrs, file=fname)
#write.table(new_ann, file = fname, sep=",",  quote=FALSE)

