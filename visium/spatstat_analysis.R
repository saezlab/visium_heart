library(Seurat)
library(SPARK)
library(tidyverse)
library(furrr)
library(FNN)

plan(multiprocess, workers = 4)
options(future.globals.maxSize= 1000*1024^2)

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

spark_analysis <- function(sample_names, tissue_paths, matrix_paths){
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
    
    write_csv(spark@res_mtest, 
              path = paste0(str_split(tissue_paths[i],"/")[[1]][1],"_spark.csv"))
    
  })
}

moran_analysis <- function(sample_names, tissue_paths, matrix_paths){
  seq(sample_names) %>% walk(function(i){
    
    d <- CreateSeuratObject(Read10X_h5(matrix_paths[i]), project = sample_names[i]) %>% SCTransform()
    
    expr <- d@assays$SCT@data %>% as.matrix %>% t %>% as.data.frame
    
    geometry <- read.csv(tissue_paths[i],
                         col.names=c("barcode","tissue","row","col","imagerow","imagecol"), 
                         header = FALSE)
    
    # just to be cautious, we will align the counts and locations ourselves 
    merged <- merge(expr %>% rownames_to_column("barcode"), 
                    geometry %>% select("barcode", "row", "col"), by = "barcode")
    
    rm(d, expr, geometry)
    
    #filter
    cat("Constructing neighborhood")
    neighbors <- knn.index(merged %>% select(row, col), k = 10) 
    neighborhood <- seq(nrow(neighbors)) %>% future_map_dfr(function(cid){
       ns <- neighbors[cid, ]
       difs <- apply(merged %>% slice(ns) %>% select(row,col), 
             1, 
             function(x) abs(x - c(merged$row[cid],merged$col[cid])))
       invalid <- unique(c(which(difs[1,]>1), which(difs[2,]>1)))
       if(!is_empty(invalid)) ns <- ns[-invalid]
       tibble(source = cid, target = ns)
    }, .progress=T)
    
    cat("Calculating Moran's I")
    moran <- seq(2, ncol(merged)-2) %>% future_map_dbl(function(i){
      y <- merged %>% pull(i)
      yavg <- mean(y)
      yvar <- var(y)
      S0 <- nrow(neighborhood)
      I <-  1/(S0*yvar) *
        (seq_along(y) %>% 
         map_dbl(~sum((y[.x] - yavg)*(y[neighborhood %>% filter(source == .x) %>% pull(target)] - yavg))) %>% 
           sum)
      I
    }, .progress = TRUE)
    
    cat("Calculating p-values")
    EI <- -1/(nrow(merged) - 1)
    
    n2  <- neighborhood
    colnames(n2) <- colnames(neighborhood)
    S1 <- 1/2 * (4*nrow(intersect(neighborhood, n2)) + nrow(setdiff(neighborhood, n2)))
    
    S2 <- nrow(merged) %>% future_map_dbl(function(i){
      nsi <- neighborhood %>% filter(source == i) %>% pull(target)
      (length(nsi) + sum(nrow(merged) %>% map_lgl(function(j){
        nsj <- neighborhood %>% filter(source == j) %>% pull(target)
        i %in% nsj
      }))^2)
    }, .progress = T) %>% sum
    
    D <- seq(2,ncol(merged)-2) %>% future_map_dbl(function(i){
      y <- merged %>% pull(i)
      yavg <- mean(y)
      sum((y - yavg)^4)/sum((y - yavg)^2)^2
    }, .progress = T)
      
    n <- nrow(merged)
    S0 <- nrow(neighborhood)
    
    A <- n*((n^2 - 3*n + 3)*S1 - n*S2 + 3*S0^2)
    B <- D*((n^2 - n)*S1 - 2*n*S2 + 6*S0^2)
    C <- (n-1)*(n-2)*(n-3)*S0^2
    
    EIsq <- (A - B)/C
    
    VI <- EIsq - EI^2
    
    zI <- (moran - EI)/sqrt(VI)
    
    pvals <- pnorm(-zI)
    adjpvals <- p.adjust(pvals, method = "fdr")
    genes <- colnames(merged)[2:(ncol(merged)-2)]
    
    result <- data.frame(gene = genes, I = moran, pvalue = pvals, adjusted_pvalue = adjpvals)
    write_csv(result, path = paste0(str_split(tissue_paths[i],"/")[[1]][1],"_moran.csv"))

  })
}

stat_aggregate <- function(ids, thr = 0.05, type = "spark"){
  ids %>% 
    map(~read_csv(paste0(type, "_results/", ., "_", type, ".csv")) %>% 
                           filter(adjusted_pvalue <= thr) %>% pull(ifelse(type=="spark",X1, gene))) %>% 
    reduce(intersect)
}

#moran_analysis(sample_names, tissue_paths, matrix_paths)

g <- read_csv("groups.csv")

hel <- stat_aggregate(g %>% pull(Healthy), type="spark")
bz <- stat_aggregate(g %>% pull(Borderzone), type="spark")
mi <- stat_aggregate(g %>% pull(MI), type="spark")
ch <- stat_aggregate(g %>% pull(Chronic), type="spark")
common <- list(hel, bz, mi, ch) %>% reduce(intersect)

result <- c("Healthy","Borderzone","MI","Chronic","Common") %>% 
  map2_dfr(list(hel,bz,mi,ch,common),~tibble(Gene = .y, Group  = .x))

result_diff <- c("Healthy","Borderzone","MI","Chronic") %>% 
  map2_dfr(list(setdiff(hel,common),setdiff(bz,common),setdiff(mi,common),setdiff(ch,common)),~tibble(Gene = .y, Group  = .x))

write_csv(result, path = "spark_results/spark_intersections.csv")
write_csv(result_diff, path = "spark_results/spark_intersections_diff.csv")



