#' In this script we manually annotate clusters to be used in subsequent
#' analyses, this happens after collapsing unneecessary states
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
              help = "tsv file with annotations and other stats related to cell filtering"),
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
  rownames_to_column("cell_id_raw") %>% 
  dplyr::mutate(cell_id = map_chr(strsplit(cell_id_raw, split = "-"), ~.x[1])) %>%
  dplyr::left_join(annotation_data, by = c("orig.ident", "cell_id"))

# Filter cells that had to be filtered

meta_data <- meta_data %>%
  dplyr::filter(!filtered_by_patspc,
                !filtered_by_undefinition)

rownames(meta_data) <- meta_data$cell_id_raw

scell_obj <- scell_obj[, meta_data$cell_id_raw]

scell_obj@meta.data <- meta_data

# Save
saveRDS(scell_obj, file = out_path)




