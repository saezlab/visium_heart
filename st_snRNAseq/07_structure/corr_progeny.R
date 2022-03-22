# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Correlate PROGENy scores with each other

library(Seurat)
library(tidyverse)

# Main

# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"

visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples) %>%
  mutate()

all_cors <- map(set_names(visium_df$visium_file, visium_df$sample), function(visium_file) { 
  print(visium_file)
  
  path_score <- readRDS(visium_file) %>%
    GetAssayData(., assay = "progeny") %>%
    as.matrix() %>%
    t() %>%
    cor()
  
  return(path_score)
})

all_cors_df <- map(all_cors, function(x) {
  
  x %>%
    as.data.frame() %>%
    rownames_to_column("feature_a") %>%
    pivot_longer(-feature_a, names_to = "feature_b")
  
  
}) %>% enframe() %>%
  unnest()

progeny_cors_plt <- all_cors_df %>%
  group_by(feature_a, feature_b) %>%
  summarise(mean_cor = mean(value)) %>%
  ggplot(., aes(x = feature_a, y = feature_b, fill = mean_cor)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
        axis.text = element_text(size = 12)) +
  xlab("") +
  ylab("") +
  coord_equal()

pdf("./results/tissue_structure/progeny_cors.pdf", height = 5, width = 5)

plot(progeny_cors_plt)

dev.off()


