# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Can I spot differential patterns of colocalization in patients?
#' 

library(tidyverse)
library(Seurat)

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


# Evidence of problems of granularity ------------------------------------------
slide_obj <- readRDS("./visium_results_manuscript/processed_visium_revisions/157781.rds")
c2l_mat <- t(as.matrix(slide_obj[["c2l_states"]]@data))
c2l <- cor(c2l_mat, method = "spearman") 
c2l <- c2l %>%
  as.data.frame() %>%
  rownames_to_column("ctA") %>%
  pivot_longer(-ctA,
               names_to = "ctB",
               values_to = "correlation")

ggplot(c2l, 
       aes(x = ctA, y = ctB, fill = correlation)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient(high = "yellow", low = "darkgreen")



# Putting cell2location scores and para/juxta representations in a single tibble ---
deconv_res <- tibble("slide_path" = slide_files_full,
                     "slide" = slide_ids) %>%
  dplyr::mutate(c2l_cors = map(slide_path, function(x) {
    
    slide_obj <- readRDS(x)
    
    c2l_mat <- t(as.matrix(slide_obj[["c2l_states"]]@data))
    
    # Filtering useless spots -----------------------
    r_ix <- rowSums(c2l_mat)
    r_ix <- r_ix >= 0.1
    c2l_mat <- c2l_mat[r_ix, ]
    
    # Correlation between cell densities ------------------------------------------
    c2l <- cor(c2l_mat, method = "spearman") 
    c2l[upper.tri(c2l, diag = T)] <- NA
    
    c2l <- c2l %>%
      as.data.frame() %>%
      rownames_to_column("ctA") %>%
      pivot_longer(-ctA,
                   names_to = "ctB",
                   values_to = "correlation") %>%
      na.omit()
  }))



deconv_res <- deconv_res %>%
  dplyr::select(-slide_path) %>%
  unnest()

deconv_res <- deconv_res %>%
  left_join(condition_dictionary, 
            by = "slide") 

top_correlations = deconv_res %>% 
  arrange(patient, -abs(correlation)) %>% 
  group_by(patient) %>%
  slice(1:50) %>% 
  dplyr::mutate(corlab = paste0(ctA, "_", ctB)) %>%
  dplyr::select(patient, corlab, correlation) %>%
  pivot_wider(names_from = patient, values_from = correlation)

deconv_res_wide <- deconv_res %>%
  dplyr::mutate(corlab = paste0(ctA, "_", ctB)) %>%
  dplyr::select(patient, corlab, correlation) %>%
  pivot_wider(names_from = patient, values_from = correlation)

deconv_res_mat <- as.matrix(deconv_res_wide[, -1])
rownames(deconv_res_mat) <- deconv_res_wide[[1]]

# Quick stats of colocalization 
top_interactions <- 0.20
variable_patterns <- sort(apply(deconv_res_mat, 1, var),decreasing = T)[1:floor(nrow(deconv_res_mat) * top_interactions)]
variable_patterns <- names(variable_patterns)

deconv_res_mat <- deconv_res_mat[variable_patterns,]

# Patient distances --------------------------------
p_dist_mat <- dist(t(deconv_res_mat))
p_hclust_data <- hclust(p_dist_mat)
p_order <- p_hclust_data$labels[p_hclust_data$order]

# Colocalization distances --------------------------------
c_dist_mat <- dist(deconv_res_mat)
c_hclust_data <- hclust(c_dist_mat)
c_order <- c_hclust_data$labels[c_hclust_data$order]
colocalization_patterns <- cutree(c_hclust_data, k = 5)
colocalization_patterns <- tibble(colocalizecells = names(colocalization_patterns),
                                  group = colocalization_patterns) %>%
  arrange(group)

# Heatmap patients --------------------------------

deconv_res %>%
  dplyr::mutate(corlab = paste0(ctA, "_", ctB)) %>%
  dplyr::filter(corlab %in% variable_patterns) %>%
  mutate(patient = factor(patient,
                          levels = p_order),
         corlab = factor(corlab,
                         levels = c_order)) %>%
  ggplot(aes(x = patient, y = corlab, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0)
  





dist_mat_plt <- as.matrix(dist_mat) %>%
  as.data.frame() %>%
  rownames_to_column("patientA") %>%
  pivot_longer(-patientA,
               values_to = "distance",
               names_to = "patientB"
              ) %>%
  dplyr::mutate(patientA = factor(patientA,
                                  levels = order),
                patientB = factor(patientB,
                                  levels = order)) %>%
  ggplot(aes(x = patientA, y = patientB, fill = distance)) +
  geom_tile()


deconv_res_mat

test = princomp((deconv_res_mat),cor = TRUE)

test = prcomp(x = t(deconv_res_mat))

test$x %>%
  as.data.frame() %>%
  rownames_to_column("patient") %>%
  ggplot(aes(x = PC1, y = PC2, label = patient)) +
  geom_point() + geom_text()

