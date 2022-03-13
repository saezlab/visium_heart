# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we explore functionally the different states/clusters
#' defined in a single seurat object

library(optparse)
library(tidyverse)
library(Seurat)
source("./utils/funcomics.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
        make_option(c("--data_path"), 
                    action ="store", 
                    default = NULL, 
                    type = 'character',
                    help = "scell seurat file"),
        make_option(c("--out_path"), 
                    action= "store", 
                    default = NULL, 
                    type = 'character',
                    help = "where to save the rds object with functional assays"),
        make_option(c("--group_class"), 
                    action= "store", 
                    default = NULL, 
                    type = 'character',
                    help = "Identity to characterize"),
        make_option(c("--gene_set_collection"), 
                    action= "store", 
                    default = NULL, 
                    type = 'character',
                    help = "path to an R object with a list of gene sets"),
        make_option(c("--module_name"), 
                    action= "store", 
                    default = NULL, 
                    type = 'character',
                    help = "name of genesets if used")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
        if(user_input=="help") next;
        cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
        assign(user_input,opt[[user_input]])
}

# Read object to analyze and define assay -----------------------------------------------------
integrated_data <- readRDS(data_path)
Idents(integrated_data) <- group_class

# Run funcomics ------------------------------------------------------------------------------

if(gene_set_collection != "not") {
        marker_dictionary <- readRDS(gene_set_collection)
} else{
        marker_dictionary <- NULL
        module_name <- NULL
}

integrated_data <- add_funcomics(visium_slide = integrated_data,
                                 species = "human",
                                 confidence_lbls = c("A","B","C","D"),
                                 top = 500,
                                 marker_dictionary = marker_dictionary,
                                 module_name = module_name,
                                 verbose = FALSE,
                                 assay = "RNA")

saveRDS(integrated_data, 
        file = out_path)
