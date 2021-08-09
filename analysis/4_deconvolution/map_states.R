#' In this script we manually annotate clusters to be used in subsequent
#' analyses
#' 
#' Warning: Hard coding is expected here, run plot_knownmarkers.R after this
#' for validation
#' 

library(optparse)
library(tidyverse)
library(Seurat)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file must have a previous annotation somewhere in meta-data, it is transformed to chr"),
  make_option(c("--dictionary_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "tsv file with annotations to transfer (only one extra column besides cell_id and orig.ident)"),
  make_option(c("--inherit_var"),
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "variable to inherit data if annotation is incomplete"),
  make_option(c("--ann_var"),
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "new variable in meta.data from dictionary_id"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds object")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Read object to annotate and get meta-data
annotation_data <- read.table(file = dictionary_path,
                              header = T,
                              sep = "\t",
                              stringsAsFactors = F)

# Read object to annotate and get meta-data, Left-join to map state annotation
scell_obj <- readRDS(data_path)
meta_data <- scell_obj@meta.data %>%
  rownames_to_column("cell_id") %>% 
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id, split = "-"), ~.x[1])) %>%
  dplyr::mutate(inherit_var = as.character(inherit_var)) %>%
  dplyr::select(orig.ident, cell_id, all_of(inherit_var)) %>%
  dplyr::left_join(annotation_data, by = c("orig.ident", "cell_id"))

meta_data[, ann_var] <- ifelse(is.na(meta_data[,ann_var]),
                               meta_data[,inherit_var],
                               meta_data[,ann_var])

# Add info to object
scell_obj <- AddMetaData(scell_obj, 
                         metadata = meta_data[, ann_var],
                         col.name = ann_var)

# Save
saveRDS(scell_obj, file = out_path)





