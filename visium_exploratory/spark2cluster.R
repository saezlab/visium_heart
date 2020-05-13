# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Mapping SPARK results to cluster annotation
#' 
#' One could query one class and extract in which clusters of which 
#' samples are the pathways overrepresented
#' 
#' 
## Function to map SPARK results to single slides
## It assumes that a given pathway that is enriched of
## SPARK genes should map back to individual slides
## in specific clusters
##
## Input:
## SPARK_annotations = result table from annotate_spark.R, run of GSA
## ora_results = tidy table of all ora results from slide_scell_map
## spark_class: class to map back to individual slides
## spark_pvalue: p_value threshold from SPARK results
## slide_pvalue: p_value threshold from ora results
## BH: use corrected p-values

SPARK2clusters = function(SPARK_annotations, ora_results, 
                          spark_class, spark_pvalue, 
                          slide_pvalue, BH = FALSE){
  
  if(spark_class == "Common"){
    sample_ids = unlist(sample_grouping)
  }else{
    sample_ids = sample_grouping[[spark_class]]
  }
  
  if(BH){
    spark_res = SPARK_annotations %>% dplyr::filter(Condition_group == spark_class,
                                                    corr_p_value <= spark_pvalue) %>%
      dplyr::select(Condition_group,gset,p_value,corr_p_value)
    
    ora_res = ora_results %>% dplyr::filter(sample_id %in% sample_ids,
                                            corr_p_value <= slide_pvalue) %>%
      dplyr::select(sample_id,gset,Cluster_ID,p_value,corr_p_value)
    
  }else{
    spark_res = SPARK_annotations %>% dplyr::filter(Condition_group == spark_class,
                                                    p_value <= spark_pvalue) %>%
      dplyr::select(Condition_group,gset,p_value,corr_p_value)
    
    ora_res = ora_results %>% dplyr::filter(sample_id %in% sample_ids,
                                            p_value <= slide_pvalue) %>%
      dplyr::select(sample_id,gset,Cluster_ID,p_value,corr_p_value)
  }
  
  all_results = left_join(spark_res, ora_res, by = "gset")
  
  colnames(all_results) = gsub(".x","_spark", colnames(all_results))
  
  colnames(all_results) = gsub(".y","_slide", colnames(all_results))
  
  return(all_results)
}

# Main
SPARK_annotations = read.csv(file = "./results/spatial_cor/spark_annotation_canonical.txt",
         sep =  "\t", header = T, stringsAsFactors = F)

sample_grouping = read.csv(file = "./results/spatial_cor/groups.csv",
                             sep =  ",", header = T, stringsAsFactors = F) %>% t()

sample_grouping = lapply(set_names(rownames(sample_grouping)), 
                         function(x){
  as.character(sample_grouping[x,])
})

# Slides with current object
slides_ids = set_names(c("157771", "157772", "157775",
               "157777", "157779","157781",
               "157782", "157785"))

ora_results = enframe(lapply(slides_ids, function(slide){
  
  dea_file = sprintf("./results/single_slide/%s/%s_differential_feat.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)$ora_canonical
  
}),name = "sample_id") %>% unnest()



spark_map = lapply(set_names(c(names(sample_grouping),"Common")), SPARK2clusters,
       SPARK_annotations = SPARK_annotations,
       ora_results = ora_results,
       spark_pvalue = 0.05,
       slide_pvalue = 0.05,
       BH = T)

# Generate objects

out_test = lapply(set_names(names(spark_map)), function(x){
  
  outfile = sprintf("./results/spatial_cor/spark_mapping_%s.tsv",
                 x)
  
  write.table(spark_map[[x]],
              file = outfile, quote = F,row.names = F,
              col.names = T, sep = "\t")
  
})






