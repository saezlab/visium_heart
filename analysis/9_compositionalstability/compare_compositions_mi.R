# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we compare the compositions of complete atlases
#' MI:
#' RNA, ATAC, Visium
#' 
#' First we show that our atlas is consistent and that per-patient, the correlation is quite high
#' 
#' 

library(tidyverse)

# Read individual compositions

atac_props <- read.csv("./results/compositions/atac_compositions.txt", sep = "\t")
rna_props <- read.csv("./results/compositions/snrna_compositions.txt", sep = "\t")
spatial_props <- read.csv("./results/compositions/spatial_compositions.txt", sep = "\t")

# Put them in a single tibble

mi_props <- left_join(spatial_props, rna_props, 
          by = c("patient_id", "cell_type")) %>%
  left_join(atac_props,
            by = c("patient_id", "cell_type")) %>%
  mutate(sn_n_cells = ifelse(is.na(sn_n_cells), 0, sn_n_cells),
         sn_prop_cells = ifelse(is.na(sn_prop_cells), 0, sn_prop_cells),
         atac_n_cells = ifelse(is.na(atac_n_cells), 0, atac_n_cells),
         atac_prop_cells = ifelse(is.na(atac_prop_cells), 0, atac_prop_cells))

# Let's correlate all combos:

mi_props <- mi_props %>%
  group_by(patient_id) %>%
  nest() %>%
  dplyr::mutate(RNAvSpatial = map(data, function(dat) {
      
      cor_res <- broom::tidy(cor.test(log1p(dat[["sp_prop_cells"]]), 
                           log1p(dat[["sn_prop_cells"]]),
                           method = "spearman"))
      
      cor_res <- cor_res[,c("estimate", "p.value")]
      colnames(cor_res) <- paste0("RNAvSpatial_",colnames(cor_res))
      cor_res
    
  })) %>%
  dplyr::mutate(RNAvATAC = map(data, function(dat) {
    
    cor_res <- broom::tidy(cor.test(log1p(dat[["atac_prop_cells"]]), 
                           log1p(dat[["sn_prop_cells"]]),
                           method = "spearman"))
    
    cor_res <- cor_res[,c("estimate", "p.value")]
    colnames(cor_res) <- paste0("RNAvATAC_",colnames(cor_res))
    cor_res
    
  })) %>%
  dplyr::mutate(SpatialvATAC = map(data, function(dat) {
    
    cor_res <- broom::tidy(cor.test(log1p(dat[["atac_prop_cells"]]), 
                           log1p(dat[["sp_prop_cells"]]),
                           method = "spearman"))
    
    cor_res <- cor_res[,c("estimate", "p.value")]
    colnames(cor_res) <- paste0("SpatialvATAC_",colnames(cor_res))
    cor_res
    
  }))


props_cors <- mi_props %>% 
  dplyr::select(-data) %>%
  unnest() %>%
  pivot_longer(-patient_id) %>%
  dplyr::mutate(comparison = strsplit(name, "_") %>%
                  map_chr(., ~ .x[[1]]),
                name = strsplit(name, "_") %>%
                  map_chr(., ~ .x[[2]])) %>%
  ungroup()


props_cors %>% dplyr::filter(name == "estimate") %>%
  group_by(comparison) %>%
  summarize(mean_cor = mean(value),
            median_cor = median(value),
            min_cor = min(value)) %>%
  write.table(., file = "./results/compositions/all_compositions_summary.txt", 
              col.names = T, row.names = F, quote = F, sep = "\t")


props_cors %>% 
  dplyr::filter(name == "estimate") %>%
  write.table(., file = "./results/compositions/all_compositions.txt", 
              col.names = T, row.names = F, quote = F, sep = "\t")

pdf("./results/compositions/all_compositions_cor.pdf", height = 3, width = 9)

props_cors %>%
  dplyr::filter(name == "estimate") %>%
  ggplot(aes(x = patient_id,
             y = comparison,
             fill = value)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom") +
  scale_fill_gradient(high = "#ffd89b", low = "#19547b", limits = c(0,1)) 

dev.off()










