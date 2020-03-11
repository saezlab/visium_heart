library(Seurat)
library(SPARK)
library(tidyverse)

sample_names <- c("AKK006 - healthy control", "AKK004 - old MI, fibrotic areas", 
                  "AKK003 - acute MI", "AKK003 - borderzone", "AKK002 - acute MI", 
                  "AKK002 - borderzone", "AKK002 - healthy part", 
                  "AKK001 - fibrotic, late stage heart failure")

tissue_paths <- c("157771/157771/outs/spatial/tissue_positions_list.txt",
                  "157772/157772/outs/spatial/tissue_positions_list.txt",
                  "157775/157775/outs/spatial/tissue_positions_list.txt",
                  "157777/157777/outs/spatial/tissue_positions_list.txt",
                  "157779/157779/outs/spatial/tissue_positions_list.txt",
                  "157781/157781/outs/spatial/tissue_positions_list.txt",
                  "157782/157782/outs/spatial/tissue_positions_list.txt",
                  "157785/157785/outs/spatial/tissue_positions_list.txt")

matrix_paths <- c("157771/157771/outs/filtered_feature_bc_matrix.h5",
                  "157772/157772/outs/filtered_feature_bc_matrix.h5",
                  "157775/157775/outs/filtered_feature_bc_matrix.h5",
                  "157777/157777/outs/filtered_feature_bc_matrix.h5",
                  "157779/157779/outs/filtered_feature_bc_matrix.h5",
                  "157781/157781/outs/filtered_feature_bc_matrix.h5",
                  "157782/157782/outs/filtered_feature_bc_matrix.h5",
                  "157785/157785/outs/filtered_feature_bc_matrix.h5")


seq(sample_names) %>% walk(function(i){

  d <- CreateSeuratObject(Read10X_h5(matrix_paths[i]), project = sample_names[i])
  counts <- d@assays$RNA@data %>% as.matrix %>% t %>% as.data.frame
  
  geometry <- read.csv(tissue_paths[i],
                       col.names=c("barcode","tissue","row","col","imagerow","imagecol"), 
                       header = FALSE)
  
  # just to be cautious, we will align the counts and locations ourselves 
  merged <- merge(counts %>% rownames_to_column("barcode"), 
                  geometry %>% select("barcode", "row", "col"), by = "barcode")
  
  spark.counts <- merged %>% select(-c("barcode","row","col")) %>% t
  colnames(spark.counts) <- merged %>% pull("barcode")
  
  #SPARK cannot handle NA's
  spark.counts[is.na(spark.counts)] <- 0
  
  spark.location <- merged %>% column_to_rownames("barcode") %>% select("row", "col")
  
  
  spark <- CreateSPARKObject(counts = spark.counts, 
                             location = spark.location,
                             percentage = 0.1, 
                             min_total_counts = 10)
    
  spark@lib_size <- apply(spark@counts, 2, sum)
  
  spark <- spark.vc(spark, 
                    covariates = NULL, 
                    lib_size = spark@lib_size, 
                    num_core = 4,
                    verbose = T)
  
  spark <- spark.test(spark, 
                      check_positive = T, 
                      verbose = T)
  
  write_csv(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")], 
            path = paste0(str_split(tissue_paths[i],"/")[[1]][1],"_spark.csv"))
  
})