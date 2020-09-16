# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here I validate TF binding activities using regulons from
#' ATAC and dorothea
#' 

library(viper)
library(tidyverse)
library(Seurat)

#' Combines two matrices in a single data frame
#' 
#' @param TF_act_matA,TF_act_matB : matrices with features in rows and ids in column 
#' @return a tibble with combined expression
#' 
#' 
combine_matrices = function(TF_act_matA, TF_act_matB){
  
  longA = pivot_longer(TF_act_matA %>% rownames_to_column("feature"),
                       -"feature",names_to = "id")
  
  longB = pivot_longer(TF_act_matB %>% rownames_to_column("feature"),
                       -"feature",names_to = "id")
  
  combined_acts = left_join(longA, longB, by = c("feature", "id")) %>%
    dplyr::rename(TFa = value.x,
                  TFb = value.y) %>% 
    na.omit()
  
  return(combined_acts)
}


#' Combines ATAC and RNA TF activities
#' Selects the top quantile from ATAC binding activities to separate
#' IDs in two groups (high, low).
#' 
#' Then I compare RNA TF activities between the two groups using Wilcox tests
#' 
#' 
#' @param TF_act_matA,TF_act_matB : matrices with features in rows and ids in columns
#' @param quantileA : quantile that estimates spots with high TF activity in A
#' @param features : features to compare if NULL the intersection
#' 
#' 
#' 

compareTFs = function(TF_act_matA, 
                      TF_act_matB, 
                      quantileA = 0.75, 
                      features = NULL,
                      plotting = F){
  
  combined_acts = combine_matrices(TF_act_matA,TF_act_matB)
  
  if(!is.null(features)){
    combined_acts = combined_acts %>%
      dplyr::filter(feature %in% features)
  }
  
  #Calculating quantiles per feature of reference set A
  
  combined_acts = combined_acts %>%
    dplyr::mutate(TF_name = feature) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(qfilt = quantile(TFa,0.75)) %>%
    dplyr::mutate(act_group =  factor(ifelse(TFa>= qfilt,
                                      "high","low"),
                                      levels = c("high","low")))
  
  wilcox_results = combined_acts %>%
    nest() %>%
    mutate(wilcox = map(data,function(df){
      stest = wilcox.test(TFb ~ act_group,data = df,
                          alternative = "g")
      broom::tidy(stest) %>%
        dplyr::mutate(corr_pvalue = p.adjust(p.value)) 
    })) %>%
    select(-data) %>% 
    unnest(wilcox) %>%
    ungroup() %>%
    dplyr::arrange(corr_pvalue)
  
    
  
  
  if(plotting){
    
    combined_acts %>%
      nest() %>%
      mutate(vplot = map(data,function(df){
        
        vplots = df %>%
          ggplot(aes(x = act_group,
                     y = TFb)) +
          geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
          theme_minimal() +
          ggtitle(unique(df$TF_name))
        
        plot(vplots)
      })) 
    
  }
 return(wilcox_results)
}

## Function to group Dorothea regulons. 
## Input: A data frame containing Dorothea regulons, as stored in 
## https://github.com/saezlab/ConservedFootprints/tree/master/data
## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}


#' Generate a regulatory network from an ATAC file
#' 
#' @param motif2gene_file: txt file with 3 columns, TF, target, confidence
#' @param cutoff: number of top genes to consider in the regulon
#' @return a viper-ready regulon
#' 

fromATACtoviper = function(motif2gene_file, cutoff = 250){
  
  regulon_df = read_delim(motif2gene_file,delim = "\t",
                          col_names = c("tf", "target", "binding_accs")) %>%
    dplyr::mutate(mor = 1, likelihood = 1) %>%
    dplyr::arrange(tf, -binding_accs) %>%
    dplyr::group_by(tf) %>%
    dplyr::slice(1:cutoff) %>%
    dplyr::ungroup() %>%
    dplyr::select(-binding_accs) %>%
    df2regulon()
  
  return(regulon_df)
  
}

getTF_matrix = function(visium_slide,viper_regulon){
  
  tf_act_mat = viper(eset = as.matrix(visium_slide[["SCT"]]@data), 
                     regulon = viper_regulon, nes = TRUE, 
                     method = "scale", minsize = 4, 
                     eset.filter = FALSE,
                     verbose = FALSE)
  
  return(tf_act_mat)
}

# Main
ATAC_to_space_validaion = function(atac, slide){
  
  if(atac %in% MotifToGene){
  
  slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                       slide,slide)
  
  # Reading the slide
  visium_slide = readRDS(slide_file)
  DefaultAssay(visium_slide) = "SCT"
  
  # HINT_TF_Activity
  hint_mat = as.matrix(visium_slide[["HINT_TF_Activity"]]@data)
  
  # Get all names
  
  atac_path = sprintf("./MotifToGene/%s",atac)
  
  cell_types = gsub("[.txt]","",list.files(atac_path)) 
  
  # Generate cell-type specific network
  allnets =  tibble("motif2gene_file" = paste0(atac_path,"/",
         list.files(atac_path)),"cell_id" = cell_types) %>%
    mutate(viper_regulon = map(motif2gene_file,fromATACtoviper)) 
  
  # Compare TF activity - takes a while
  allnets = allnets %>%
    mutate(comp_results = map(allnets$viper_regulon, function(vp){
      
      TF_act_matB = getTF_matrix(visium_slide,vp)
      
      return(compareTFs(TF_act_matA = as.data.frame(hint_mat),
                        TF_act_matB = as.data.frame(TF_act_matB)))
    }))
  
  saveRDS(allnets, 
          file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac))
  
  return(allnets)
  
  }
  
}

# Healthy

sample_dictionary = read_delim("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                               delim = "\t") %>%
  dplyr::select(snATAC,Visium) %>%
  mutate_all(as.character)

MotifToGene = list.files("./MotifToGene/")

TFcomparison = sample_dictionary[,] %>% 
  mutate(TF_comp = map2(snATAC,Visium,ATAC_to_space_validaion))


# Specific analysis throughout the paper


# Healthy - cardio 
# CK166

atac = "CK166"
slide = "157771"

# Read visium
slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                     slide,slide)
visium_slide = readRDS(slide_file)
DefaultAssay(visium_slide) = "SCT"

# Getting the regulatory networks
atac_path = sprintf("./MotifToGene/%s",atac)
all_res = readRDS(file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac))

cardio_TFs = all_res$comp_results[[1]] %>%
  dplyr::filter(corr_pvalue < 0.01) %>%
  dplyr::slice(1:10) %>% pull(feature)

TF_act_matB = getTF_matrix(visium_slide,all_res$viper_regulon[[1]][cardio_TFs])
hint_mat = as.matrix(visium_slide[["HINT_TF_Activity"]]@data)
dorothea_mat = as.matrix(visium_slide[["dorothea"]]@data)
  
compareTFs(TF_act_matA = as.data.frame(hint_mat),
           TF_act_matB = as.data.frame(TF_act_matB),
           features = cardio_TFs,plotting = T)

compareTFs(TF_act_matA = as.data.frame(hint_mat),
           TF_act_matB = as.data.frame(dorothea_mat),
           features = cardio_TFs,plotting = T)



visium_slide[['cardioTFs']] = CreateAssayObject(data = TF_act_matB)
DefaultAssay(visium_slide) = "cardioTFs"
SpatialFeaturePlot(visium_slide,features = "MEF2C")


allnets = allnets %>%
  mutate(comp_results = map(allnets$viper_regulon, function(vp){
    
    TF_act_matB = getTF_matrix(visium_slide,vp)
    
    return(compareTFs(TF_act_matA = as.data.frame(hint_mat),
                      TF_act_matB = as.data.frame(TF_act_matB)))
  }))



sample_dictionary













