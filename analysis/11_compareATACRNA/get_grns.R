# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we build networks and evaluate them before comparing them with HINT results
#' All networks have the same number of TFs but not interactions

grns <- list.files("./reg_nets/raw/")

grns <- set_names(paste0("./reg_nets/raw/", grns), 
                  gsub("[.]txt","",grns)) 

all_nets <- map(grns, read_table2, col_names = F) %>%
  enframe() %>%
  unnest()

colnames(all_nets) <- c("cell_type", "source", 
                        "target", "n_reads")

# We will weight them by the max number of reads around the binding site

w_grns <- all_nets %>%
  group_by(cell_type, source) %>%
  mutate(max_reads = max(n_reads)) %>%
  mutate(likelihood = n_reads/max_reads) %>%
  dplyr::filter(likelihood > 0.3) %>%
  dplyr::mutate(mor = 1) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(out_file = paste0("./reg_nets/processed/", cell_type,".txt"))

walk2(w_grns$data, w_grns$out_file, function(dat, f){
  print(f)
  write.table(dat, file = f,col.names = T, row.names = F, sep = "\t", quote = F)
})

