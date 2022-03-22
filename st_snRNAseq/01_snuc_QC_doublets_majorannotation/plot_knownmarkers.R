# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we automatize visualization of commonly known heart markers
#' 
#' It requires a folder path with processed Seurat objects from run_singleprocessing.R

library(optparse)
library(tidyverse)
library(Seurat)
source("./analysis/utils/sc_plts.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--folder"), 
              action = "store_true", 
              default = FALSE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%sample.rds"),
  make_option(c("--downsampling"), 
              action = "store_true", 
              default = FALSE, 
              type = 'logical',
              help = "should the object subsampled to the 50%?"), 
  make_option(c("--used_assay"), 
              action ="store", 
              default = "RNA", 
              type = 'character',
              help = "name of the sample if single path provided"),
  make_option(c("--id_label"), 
              action ="store", 
              default = "opt_clust", 
              type = 'character',
              help = "name of the sample if single path provided"),
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "direct folder with .rds files or .rds files"),
  make_option(c("--out_fig_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "directory where to save the pdf file")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# This is an option to give a folder and parse all the files to process --------------------------------
if(folder) {
  slide_files <- list.files(path)
  slide_files_path <- paste0(path,slide_files)
} else {
  slide_files <- unlist(strsplit(path,"/")) %>% 
    dplyr::last()
  slide_files_path <- path
}

out_file <- gsub(pattern = "rds",
                 replacement = "pdf",
                 slide_files)

# Make data frame of parameters ------------------------------------------------------------------------
param_df <- tibble(slide_file = slide_files_path,
                   out_fig_file = paste0(out_fig_path, "markers_", out_file))

# Read marker list ------------------------------------------------------------------------------------

# Gene markers from Christoph
markers_stable <- (read.table("/beegfs/work/hd_wh241/MI_revisions/markers/Kuppe_def.txt",
                              sep = "\t",stringsAsFactors = F,
                              header = T))[,c(1,2)]

markers_stable[,1] <- toupper(markers_stable[,1])
colnames(markers_stable) <- c("gene","cell_type")

# Gene markers from Hubner

alt_markers <- list(
  fibroblasts = c("DCN","GSN","PDGFRA"),
  pericytes = c("RGS5","ABCC9","KCNJ8"),
  adipocytes = c("GPAM","LEP"),
  smooth_muscle = c("MYH11", "TAGLN", "ACTA2"),
  neuronal = c("PLP1","NRXN1","NRXN3"),
  endothelial = c("VWF","PECAM1","CDH5"),
  atrial_cardio = c("NPPA","MYL7","MYL4"),
  ventricular_cardio = c("MYH7","MYL2","PRDM16"),
  myeloid = c("CD14","C1QA","FOLR2"),
  lymphoid = c("CD8A","NKG7","LCK")
) %>% enframe(name = "cell_type",
              value = "gene") %>%
  unnest() %>%
  dplyr::select(gene, cell_type) %>%
  as.data.frame()

# Gene markers from Christoph but coming from first iteration of data

ck_markers <- read_table2("/beegfs/work/hd_wh241/MI_revisions/markers/CK_mi_markers.txt") %>%
  as.data.frame()

# Make wrapper function

plt_markers <- function(slide_file, out_fig_file) {
  
  scell_obj <- readRDS(slide_file)
  
  if(downsampling == TRUE) {
    
    print("downsampling to 50%")
    
    set.seed(241099)
    
    cells <- scell_obj@meta.data %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(cell_id, id_label) %>%
      group_by_at(id_label) %>%
      nest() %>%
      mutate(data = map(data, ~.x[[1]])) %>%
      mutate(n_cells = map(data, length)) %>%
      unnest(n_cells) %>%
      mutate(n_cells = floor(n_cells * 0.5)) %>%
      mutate(sampled_cells = map2(data, n_cells, function(x, y) {
        sample(x, y)
      })) %>%
      pull(sampled_cells) %>%
      unlist()
    
    scell_obj <- scell_obj[, cells]
    
  }
  
  
  home_hmap <- domarker_hmap(SeuratObject = scell_obj,
                assay = used_assay,
                identity_label = id_label, 
                GSC = markers_stable)
  
  hubner_hmap <- domarker_hmap(SeuratObject = scell_obj,
                             assay = used_assay,
                             identity_label = id_label, 
                             GSC = alt_markers)
  
  ck_hmap <- domarker_hmap(SeuratObject = scell_obj,
                               assay = used_assay,
                               identity_label = id_label, 
                               GSC = ck_markers)
  
  pdf(out_fig_file, width = 19, height = 15)
  print(home_hmap)
  print(hubner_hmap)
  print(ck_hmap)
  dev.off()
  
}

# Make data frame of parameters ------------------------------------------------------------------------
pwalk(param_df, plt_markers)





