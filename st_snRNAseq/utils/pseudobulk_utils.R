# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Utilities to analyze single cell pseudobulked data
#' 
library(edgeR)

#' Filter expression matrix
#' @param expression_matrix: count matrix with samples in columns and genes in rows
#' @param min.count: numeric. Minimum count required for at least some samples.
#' @param min.prop: numeric. Minimum proportion of samples in the smallest group that express the gene.

edgeR_filtering <- function(expression_matrix, 
                            min.count = 10,
                            min.prop = 0.1,
                            min.total.count = 15) {
  
  meta_data <- data.frame("sample_names" = colnames(expression_matrix),
                          "lib.size" = colSums(expression_matrix),
                          stringsAsFactors = F)
  
  #edgeR pipeline
  bulk_data <- edgeR::DGEList(expression_matrix, 
                      samples= meta_data)
  
  #Filtering lowly expressed genes in the context of the groups
  keep = edgeR::filterByExpr(bulk_data, 
                             min.count = min.count,
                             min.prop = min.prop,
                             min.total.count = min.total.count,
                             group = rep("same_group",ncol(bulk_data)))
  
  bulk_data = bulk_data[keep,]
  
  
  return(bulk_data$counts)
}

#' Counts per million
#' @param expression_matrix: count matrix with samples in columns and genes in rows
#' @param scale_factor: used to generate the final counts
cpm_norm <- function(expression_matrix, scale_factor = 10000) {
  # First we transpose because of R conventions
  # Correct by coverage
  norm_mtx <- t(t(expression_matrix)/colSums(expression_matrix))
  norm_mtx <- log1p(norm_mtx * scale_factor)
  return(norm_mtx)
}


# Implementation of distances
#' @param matA: count matrix with samples in columns and genes in rows
#' @param matB: count matrix with samples in columns and genes in rows
#' @param jsd: use Jensen-Shannon divergence distance? (only for raw counts)
compare_profiles <- function(matA, matB, jsd = TRUE, matA_label = "A", matB_label = "B") {
  
  gene_ids <- intersect(rownames(matA), 
                        rownames(matB))
  
  colnames(matA) <- paste0(matA_label, colnames(matA))
  
  colnames(matB) <- paste0(matB_label, colnames(matB))
  
  matA <- matA[gene_ids, colSums(matA)>0]
  matB <- matB[gene_ids, colSums(matB)>0]
  
  if(jsd){
    
    dist_mat <- philentropy::JSD(t(cbind(matA, matB)), est.prob = "empirical")
    rownames(dist_mat) = colnames(dist_mat) <- c(colnames(matA),colnames(matB))
    dist_mat <- dist_mat ** (1/2)
    dist_mat <- dist_mat[colnames(matA),colnames(matB)]
    
    return(dist_mat)
    
  } else {
    
    dist_mat <- cor((cbind(matA, matB)), method = "spearman")
    dist_mat <- dist_mat[colnames(matA),colnames(matB)]
    
    return(dist_mat)
  }
}

# separation of matrices per cell-type
#' @param mat: Pseudobulk expression matrix
#' @param cell_state_dic: Dictionary to map cell-type to cell_states

separate_pseudobulk = function(mat, cell_state_dic) { 
  
  colnames(cell_state_dic) <- c("major", "cell_states")
  
  long_mat <- as.data.frame(mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, values_to = "counts",
                 names_to = "cell_states") %>%
    left_join(cell_state_dic) %>%
    group_by(major) %>%
    nest() %>%
    mutate(expr_mat = map(data, function(x) {
      
      temp_mat <- pivot_wider(x,
                              names_from = cell_states,
                              values_from = counts)
      
      rows_ids <- temp_mat[[1]]
      temp_mat <- as.matrix(temp_mat[, -1])
      rownames(temp_mat) <- rows_ids
      return(temp_mat)
      
    }))
  
  return(long_mat)
  
}


#' Filter matrix by gene
#' @param expression_mat: a expression matrix with genes in rows and samples in colums
#' @param gene_list: do you want to use jensen shannon divergences distances?
#' @return a reduced matrix containing only genes of interest
filter_genes <- function(expression_mat, gene_list) {
  
  gene_list <- gene_list[gene_list %in% rownames(expression_mat)]
  
  if(length(gene_list) > 0) {
    
    return(expression_mat[gene_list,])
    
  } else{
    
    NULL
    
  }
}

#' Compare expression profiles of the same matrix
#' @param expression_mat: a expression matrix with genes in rows and samples in colums
#' @param jsd: do you want to use jensen shannon divergences distances?
#' @return a distance matrix
self_compare <- function(expression_mat, jsd = T) {
  
  if(jsd) {
    compare_profiles(expression_mat,
                     expression_mat)
    
  } else {
    
    compare_profiles(expression_mat,
                     expression_mat,
                     jsd = F)
    
  }
  
}












