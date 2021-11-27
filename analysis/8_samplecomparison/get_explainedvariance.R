# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Calculate the amount of explained variance that can be assigned to 
#' the phenotypic differences in each cell-type

library(scater)
library(tidyverse)
library(rdist)
source("./analysis/utils/pseudobulk_utils.R")

all_mats <- readRDS("./processed_snrnaseq/integration/psxpat_integrated_rnasamples_filt.rds")
pat_anns <- readRDS("./markers/snrna_patient_anns_revisions.rds")[, c("patient_id", "patient_group")] %>%
  unique()

get_propexplvar <- function(count_matrix, pat_anns) {
  
  norm_mat <- cpm_norm(count_matrix)
  pcs <- prcomp(x = t(norm_mat)) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_id") %>%
    left_join(pat_anns)
  
  plot(pc_summary %>%
    ggplot(aes(x = PC1, y = PC2, color = patient_group)) +
    geom_point())
  
  pc_summary <- pc_summary %>%
    pivot_longer(-c("patient_id", "patient_group")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(aov_res = map(data, function(dat) {
      aov(value ~ patient_group,data = dat) %>%
        broom::tidy()
    }))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x))
  
  pc_summary <- pc_summary %>% 
    dplyr::select(aov_res) %>%
    unnest() %>%
    dplyr::filter(term != "Residuals") %>%
    left_join(prop_var) %>%
    ungroup() %>%
    dplyr::mutate(p.adj = p.adjust(p.value))
    
}

expl_vars <- all_mats %>%
  dplyr::mutate(aov_res = map(count_matrix, get_propexplvar, pat_anns = pat_anns))


expl_var_df <- expl_vars %>%
  dplyr::select(cell_type, aov_res) %>%
  unnest() %>%
  dplyr::mutate(p.adj = p.adjust(p.value)) %>%
  dplyr::filter(p.adj < 0.15) %>%
  group_by(cell_type) %>%
  summarise(expl_var = sum(expl_var)) %>%
  arrange(-expl_var)

ct_order <- pull(expl_var_df, cell_type)

expl_var_plt <- ggplot(expl_var_df, aes(x = factor(cell_type,
                                   levels = ct_order),
                        y = expl_var)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust =0.5),
        axis.title.y = element_text(size = 12)) +
  ylab("Prop. explained variance") +
  xlab("")
  
write_csv(expl_var_df, file = "./results/sample_comparison/expl_var_GEX.csv")

pdf("./results/sample_comparison/expl_var_GEX.pdf", height = 3.5, width = 3.5)

plot(expl_var_plt)

dev.off()




