# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we build a tissue consensus network
#' 

grns <- list.files("./reg_nets/processed/")
grns <- set_names(paste0("./reg_nets/processed/", grns), 
                  gsub("[.]txt","",grns)) 

all_nets <- map(grns, read_table2, col_names = T) %>%
  enframe() %>%
  unnest()

norm_factor <- unique(all_nets$name) %>%
  length()

all_nets <- all_nets %>%
  mutate() %>%
  dplyr::select(-c("name", "n_reads", "max_reads")) %>%
  group_by(source, target) %>%
  dplyr::summarize(likelihood = sum(likelihood)) %>%
  dplyr::mutate(likelihood = likelihood/norm_factor) %>%
  dplyr::mutate(mor = 1)

saveRDS(all_nets, file = "./reg_nets/tissue_consensus_net.rds")










