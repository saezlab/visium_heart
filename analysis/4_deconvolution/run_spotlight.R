# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Deconvolute slides using spotlight and a set of DEGs from cell-states or
#' any grouping variable
library(optparse)
library(Seurat)
library(tidyverse)
library(SPOTlight)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--scell_data"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell data with states in a variable"),
  make_option(c("--ident_var"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "variable to group data"),
  make_option(c("--folder"), 
              action = "store_true", 
              default = TRUE, 
              type = 'logical',
              help = "is the path added a folder with structure ./%slides.rds"),
  make_option(c("--data_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "path with visium data (r objects)"),
  make_option(c("--out_fig_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the pdf objects with qc?"),
  make_option(c("--out_file"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the deconvolution matrix?")

)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# This is an option to give a folder and parse all the files to process --------------------------------
if(folder) {
  sample_names <- list.files(data_path)
  slide_files <- paste0(data_path,sample_names)
  sample_names <- gsub(pattern = "[.]rds",
                       replacement = "",
                       sample_names)
} else {
  slide_files <- data_path
  sample_names <- gsub(pattern = "[.]rds",
                       replacement = "",
                       strsplit(data_path, "/") %>% last())
}

# Create a tibble with run parameters --------------------------------

param_df <- tibble(sample_name = sample_names,
                   slide_file = slide_files,
                   out_file = paste0(out_file, sample_names, "_deconvmat.rds"),
                   out_fig_file = paste0(out_fig_file, sample_names, "_spotlightqc.pdf"))

# Spotlight pipeline --------------------------------

# Processing scell data -----------------------------------------------------
scell_dat <- readRDS(file = scell_data)
scell_dat$deconv_state <- paste0("state_",
                                 scell_dat@meta.data[,ident_var])
Idents(scell_dat) <-  "deconv_state"
ct_markers <- FindAllMarkers(scell_dat,
                             only.pos = TRUE, 
                             logfc.threshold = 0.5,
                             min.pct = 0.5)

# Filter data so that only cells with n-differentially expressed genes can be deconvoluted
n = 10

states <- ct_markers %>%
  group_by(cluster) %>%
  summarize(ngenes = length(cluster)) %>%
  dplyr::filter(ngenes >= n) %>%
  pull(cluster) %>% as.character()

ct_markers <- dplyr::filter(ct_markers,
                            cluster %in% states)

scell_dat <- subset(scell_dat, 
                     subset = deconv_state %in% states)

run_spot_pplne <- function(sample_name, slide_file, out_file, out_fig_file) {
  
  print(sample_name)
  
  set.seed(123)
  
  # Running spotlight ~ 10 mins -----------------------------------------------------
  visium_slide <- readRDS(file = slide_file)
  
  spotlight_ls <- spotlight_deconvolution(se_sc = scell_dat,
                                          counts_spatial = visium_slide@assays$Spatial@counts,
                                          clust_vr = "deconv_state",
                                          cluster_markers = ct_markers,
                                          method = "nsNMF",
                                          min_cont = 0.05)
  
  print("deconv completed")
  
  decon_mtrx <- spotlight_ls[[2]]
  decon_mtrx <- decon_mtrx[,!grepl("res_ss", colnames(decon_mtrx))]
  rownames(decon_mtrx) <- rownames(visium_slide@meta.data)
  colnames(decon_mtrx) <- paste0(gsub("[.]", "-", colnames(decon_mtrx)))
  
  nmf_mod <- spotlight_ls[[1]]
  h <- NMF::coef(nmf_mod[[1]])
  rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
  
  topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
    h = h,
    train_cell_clust = nmf_mod[[2]])
  
  topics_profiles <- topic_profile_plts[[2]] + 
    theme(axis.text.x = element_text(angle = 90), 
          axis.text = element_text(size = 12))
  
  pdf(out_fig_file, height = 12, width = 13)
  
  print(topics_profiles)
  
  dev.off()
  
  saveRDS(object = decon_mtrx, file = out_file)
  
}

# Main

pwalk(param_df, run_spot_pplne)











