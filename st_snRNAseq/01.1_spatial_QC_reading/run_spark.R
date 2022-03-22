# Copyright (c) [2021] [Jovan Tanevski, Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we run a SPARK analysis from a folder of the following structure:
#' folder
#' |
#' sample.rds

library(optparse)
library(tidyverse)
library(Seurat)
library(SPARK)
library(furrr)

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "direct filtered_feature_bc_matrix file"),
  make_option(c("--out_path"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "where to save the rds objects")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list=option_list))

cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Make data frame of parameters ------------------------------------------------------------------------
param_df <- tibble(sample_names = list.files(path) %>%
                     gsub(".rds", "", .),
                   tissue_paths =  list.files(path,full.names = T)) %>%
  mutate(out_file = paste0(out_path, sample_names, "_spark.csv"))

# Automatic spark processing of each data set --------------------------------

spark_analysis <- function(sample_names, tissue_paths, out_file){
    print(sample_names)
    sample_seurat <- readRDS(tissue_paths)
    
    counts <- GetAssayData(sample_seurat, assay = "Spatial") %>% as.matrix %>% t %>% as.data.frame
    
    geometry <- GetTissueCoordinates(sample_seurat,
                                    cols = c("row", "col"),
                                    scale = NULL) %>%
      rownames_to_column("barcode")
    
    # just to be cautious, we will align the counts and locations ourselves 
    merged <- merge(counts %>% rownames_to_column("barcode"), 
                    geometry %>% select("barcode", "row", "col"), by = "barcode")
    
    # delete unnecesary objects
    rm(sample_seurat)
    rm(counts)
    rm(geometry)
    
    spark.counts <- merged %>% select(-c("barcode","row","col")) %>% t
    colnames(spark.counts) <- merged %>% pull("barcode")
    
    #SPARK cannot handle NA's
    spark.counts[is.na(spark.counts)] <- 0
    
    genes_nreads_total <- (rowSums(spark.counts > 0)/ncol(spark.counts))
    genes_nreads_total <- genes_nreads_total[genes_nreads_total >= 0.01]
    
    spark.counts <- spark.counts[names(genes_nreads_total), ]
    
    #Let's reduce the feature space (genes with super low counts)
    spark.location <- merged %>% column_to_rownames("barcode") %>% select("row", "col")
    
    rm(merged)
    
    sparkX <- sparkx(spark.counts,
                     spark.location,
                     numCores=4,
                     option="mixture")
    
    write_csv(sparkX$res_mtest %>% 
                rownames_to_column("gene"), 
              path = out_file, col_names = T)

}

# Main -----------------

pwalk(param_df, spark_analysis)

