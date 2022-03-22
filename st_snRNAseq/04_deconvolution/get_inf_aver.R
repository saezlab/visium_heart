# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' From cell2location estimates, get a summarized inferred average
#' 

library(tidyverse)

nb_estimates_folder <- "./results/nb_estimates_csv/"
nb_files <- set_names(paste0(nb_estimates_folder ,
                             list.files(nb_estimates_folder)), 
                      paste0("iter_", seq(1,5,1)))

nb_mats <- map(nb_files, function(x) {
  read_csv(x, col_names = T) %>%
    dplyr::rename("gene" = X1)
})


nb_mats <- nb_mats %>%
  enframe() %>%
  unnest() %>%
  pivot_longer(-c(name, gene), 
               names_to = "cell_type", 
               values_to = "infer_aver")


map(set_names(nb_mats$cell_type %>% unique), function(ct) {
  
  ct_dat <- nb_mats %>%
    dplyr::filter(cell_type == ct) %>%
    dplyr::select(-cell_type) %>%
    pivot_wider(values_from = infer_aver,
                names_from = name,
                values_fill = NA) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    na.omit()
  
  cor_ct_dat <- cor(ct_dat)
  
  mean(cor_ct_dat[upper.tri(cor_ct_dat,diag = F)])
  
}) %>% enframe() %>%
  unnest() %>%
  ggplot(aes(x = value, y = name)) +
  geom_bar(stat = "identity")


# Final data 

nb_mats <- nb_mats %>%
  mutate(cell_type = gsub("-", "_", cell_type)) %>%
  group_by(gene, cell_type) %>%
  summarize(mean_infer_aver = mean(infer_aver)) %>%
  ungroup() %>%
  pivot_wider(names_from = cell_type, values_from = mean_infer_aver) %>%
  column_to_rownames("gene")

write.csv(nb_mats, file = paste0(nb_estimates_folder, "mean_infer_aver.csv"), row.names = T, col.names = T, quote = F)









