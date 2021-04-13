# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' How do cells correlate with functional scores?

library(tidyverse)
library(Seurat)
library(CCA)
library(CCP)

# Get all data files ------------------------------------------------------------
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)
slide_files_full <- paste0(slide_files_folder, 
                           slide_files)
slide_ids <- gsub("[.]rds", "", slide_files)

# Annotate slides -------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(Visium, New.Ids) %>%
  dplyr::rename(slide = Visium,
                patient = New.Ids)


# Putting cell2location scores and para/juxta representations in a single tibble ---
deconv_res <- tibble("slide_path" = slide_files_full,
                     "slide" = slide_ids) %>%
  dplyr::mutate(assay_mats = map(slide_path, function(x) {
    
    slide_obj <- readRDS(x)
    
    id_order <- colnames(slide_obj)
    
    assay_mats <- list("c2l" = t(as.matrix(slide_obj[["c2l"]]@data))[id_order,],
                       "progeny" = t(as.matrix(slide_obj[["progeny"]]@data))[id_order,])
  }))

deconv_res <- deconv_res %>%
  left_join(condition_dictionary, 
            by = "slide")

# Performing univariate correlations ----------------------------
univariate_cor <- deconv_res %>%
  dplyr::select(slide, assay_mats) %>%
  dplyr::mutate(univariate_cor = map(assay_mats, function(x) {
    
    cor_res <- cor(x[["c2l"]], x[["progeny"]]) %>%
      as.data.frame() %>%
      rownames_to_column("cell_type") %>%
      pivot_longer(-cell_type,
                   values_to = "correlation",
                   names_to = "pathway")
  })) %>%
  dplyr::select(slide, univariate_cor, patient) %>%
  unnest() %>%
  dplyr::mutate(corlabel = paste0(cell_type, "_",pathway))

cor_order <- univariate_cor %>%
  group_by(corlabel) %>%
  summarise(mediancor = median(correlation)) %>%
  arrange(-abs(mediancor)) %>%
  slice(1:30) %>%
  arrange(mediancor) %>%
  pull(corlabel)

univariate_plt <- univariate_cor %>%
  dplyr::filter(corlabel %in% cor_order) %>%
  mutate(corlabel = factor(corlabel,
                           levels = cor_order)) %>%
  ggplot(aes(y = corlabel, 
             x = correlation)) +
  geom_boxplot() +
  geom_point(aes(color = patient))
  
pdf("visium_results_manuscript/structure/PROGENy_c2l_univariatecors.pdf", 
    height = 8, width = 7)

plot(univariate_plt)

dev.off()

# Performing multivariate correlations ----------------------------


#' Perform CCA: correlates columns of matrices
#' @param mata: data matrix with continous values
#' @param matb: data matrix with continous values
#' @return a data frame with standarized coefficients for significant canonical dimensions
cca_cor <- function(mata, matb) {
  
  # Standardize columns
  mata <- scale(mata, center = T, scale = T)
  matb <- scale(matb, center = T, scale = T)
  
  # Run CCA
  cc1 <- cc(mata, matb)
  cc2 <- comput(mata, matb, cc1)
  
  # tests of canonical dimensions
  rho <- cc1$cor
  # Define number of observations, number of variables in first set, and number of variables in the second set.
  n <- dim(mata)[1]
  p <- ncol(mata)
  q <- ncol(matb)
  
  # Calculate p-values using the F-approximations of different test statistics:
  sign_test <- p.asym(rho, n, p, q, tstat = "Wilks")
  n_canonical <- (sign_test$p.value < 0.05)
  
  # Standardized coefficients
  s1 <- diag(sqrt(diag(cov(mata))))
  mata_scoefs <- s1 %*% cc1$xcoef
  rownames(mata_scoefs) <- colnames(mata)
  colnames(mata_scoefs) <- paste0("CD",1:ncol(mata_scoefs))
  mata_scoefs <- mata_scoefs[, n_canonical]
    
  s2 <- diag(sqrt(diag(cov(matb))))
  matb_scoefs <- s2 %*% cc1$ycoef
  rownames(matb_scoefs) <- colnames(matb)
  colnames(matb_scoefs) <- paste0("CD",1:ncol(matb_scoefs))
  matb_scoefs <- matb_scoefs[, n_canonical]
  
  # Final coefficients
  
  hclust_res <- hclust(dist(rbind(mata_scoefs, 
               matb_scoefs)))
  
  order_feats <- hclust_res$labels[hclust_res$order]
  
  all_res <- rbind(mata_scoefs, 
                   matb_scoefs) %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, 
                 values_to = "std_coef",
                 names_to = "canonical_dimension") %>%
    mutate(feature = factor(feature,
                            levels = order_feats))
  
  return(list("can_corr_coef" = rho[n_canonical],
              "std_coefs" = all_res))
}

# Using all samples:

all_progeny <- map(deconv_res$assay_mats, ~.x[["progeny"]])
all_progeny <- reduce(all_progeny, rbind)

all_c2l <- map(deconv_res$assay_mats, ~.x[["c2l"]])
all_c2l <- reduce(all_c2l, rbind)

cca_all <- cca_cor(all_c2l, all_progeny)

pdf("./visium_results_manuscript/structure/cca_progeny_c2l.pdf", height = 6, width = 6)

print(ggplot(cca_all$std_coefs, aes(x = canonical_dimension,
                    y = feature,
                    fill = std_coef)) +
  geom_tile() +
  scale_fill_gradient2())

dev.off()

# Filtering low quality:
deconv_res <- deconv_res %>%
  dplyr::filter(! patient %in% c("P2_IZ", "P3_IZ"))

all_progeny <- map(deconv_res$assay_mats, ~.x[["progeny"]])
all_progeny <- reduce(all_progeny, rbind)

all_c2l <- map(deconv_res$assay_mats, ~.x[["c2l"]])
all_c2l <- reduce(all_c2l, rbind)

cca_all <- cca_cor(all_c2l, all_progeny)

pdf("./visium_results_manuscript/structure/cca_progeny_c2l_woiz.pdf", height = 6, width = 6)

ggplot(cca_all$std_coefs, aes(x = canonical_dimension,
                              y = feature,
                              fill = std_coef)) +
  geom_tile() +
  scale_fill_gradient2()

dev.off()
