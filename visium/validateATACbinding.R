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
                      plotting = F,
                      path = NULL){
  
  combined_acts = combine_matrices(TF_act_matA,TF_act_matB)
  
  if(!is.null(features)){
    combined_acts = combined_acts %>%
      dplyr::filter(feature %in% features)
  }
  
  #Calculating quantiles per feature of reference set A
  
  combined_acts = combined_acts %>%
    dplyr::mutate(TF_name = feature) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(qfilt = quantile(TFa,quantileA)) %>%
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
    plot_list = combined_acts %>%
      nest() %>%
      mutate(vplot = map(data,function(df){
        
        vplots = df %>%
          ggplot(aes(x = act_group,
                     y = TFb)) +
          geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
          theme_minimal() +
          ggtitle(unique(df$TF_name))
      })) %>%
      dplyr::pull(vplot)
    
    pdf(file = path, width = 5,height = 8)
    
    map(plot_list,plot)
    
    dev.off()
    
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

#' Use viper regulatory network and make it ready for Module Score
#' 
#' @param motif2gene_file: txt file with 3 columns, TF, target, confidence
#' @param cutoff: number of top genes to consider in the regulon
#' @return a viper-ready regulon
#' 

fromvipertoMS = function(viper_regulon){
  
  MS_regulon = map(viper_regulon, function(x){
    return(names(x$tfmode))
  })
  
  return(MS_regulon)
}


#' Generate a regulatory network from an ATAC file: TO DO
#' 
#' @param motif2gene_file: txt file with 3 columns, TF, target, confidence
#' @param cutoff: number of top genes to consider in the regulon
#' @return a viper-ready regulon
#' 
getTF_matrix = function(visium_slide,viper_regulon){
  
  tf_act_mat = viper(eset = as.matrix(visium_slide[["SCT"]]@data), 
                     regulon = viper_regulon, nes = TRUE, 
                     method = "scale", minsize = 4, 
                     eset.filter = FALSE,
                     verbose = FALSE)
  
  return(tf_act_mat)
}


#' Generate a regulatory network from an ATAC file: TO DO
#' 
#' @param motif2gene_file: txt file with 3 columns, TF, target, confidence
#' @param cutoff: number of top genes to consider in the regulon
#' @return a viper-ready regulon
#' 
getTF_matrix_MS = function(visium_slide,MS_regulon){
  
  names_vect = gsub("[.]","_",names(MS_regulon))
  names_vect = gsub("-","_",names_vect)
  
  tf_act_mat = AddModuleScore(visium_slide,features = MS_regulon,
                              name = paste0(names_vect,"__"))
  
  tf_act_mat = tf_act_mat@meta.data
  
  cell_ids = rownames(tf_act_mat)
  calculated_regulons = colnames(tf_act_mat)[grepl("__", colnames(tf_act_mat))]
  
  tf_act_mat = tf_act_mat[,calculated_regulons]
  
  colnames(tf_act_mat) = unlist(map(strsplit(colnames(tf_act_mat),split = "__"), function(x) x[1]))
  
  rownames(tf_act_mat) = cell_ids
  
  tf_act_mat = t(as.matrix(tf_act_mat))
  
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
    mutate(viper_regulon = map(motif2gene_file, fromATACtoviper)) %>%
    mutate(MS_regulon = map(viper_regulon, fromvipertoMS))
  
  # Compare TF activity - takes a while
  allnets = allnets %>%
    mutate(comp_results = map(allnets$MS_regulon, function(vp){
      
      TF_act_matB = getTF_matrix_MS(visium_slide,vp)
      
      return(compareTFs(TF_act_matA = as.data.frame(hint_mat),
                        TF_act_matB = as.data.frame(TF_act_matB)))
    }))
  
  saveRDS(allnets, 
          file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac))
  
  return(allnets)
  
  }
  
}

# All runs

sample_dictionary = read_delim("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                               delim = "\t") %>%
  dplyr::select(snATAC,Visium) %>%
  mutate_all(as.character)

MotifToGene = list.files("./MotifToGene/")

TFcomparison = sample_dictionary %>% 
  mutate(TF_comp = map2(snATAC,Visium,ATAC_to_space_validaion))

# Visualizations

# Healthy
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
all_res = readRDS(file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac)) %>%
  mutate(MS_regulon = map(MS_regulon, function(x){
    names(x) = toupper(names(x))
    return(x)
  }))

CM_TFs = c("MEF2C", "NR5A2")
ENDO1_TFs = c("SOX8")
FIB1_TFs = c("TCF21")
FIB2_TFs = c("FOS", "JUN")
MACROPH_TFs = c("ETV6","SPIB")

TFs = c(CM_TFs, ENDO1_TFs, FIB1_TFs, FIB2_TFs, MACROPH_TFs)

hint_mat = as.matrix(visium_slide[["HINT_TF_Activity"]]@data)

all_res = all_res %>%
  dplyr::filter(cell_id %in% c("Cardiomyocyes",
                               "Endohelial_1",
                               "Fibroblass_1",
                               "Fibroblass_2",
                               "Macrophages_1")) %>%
  group_by(cell_id) %>%
  dplyr::mutate(TF_act_matB = map(MS_regulon, function(x){
    
    TF_act = getTF_matrix_MS(visium_slide,MS_regulon = x[TFs])
    
    return(TF_act)
  }))

all_res$TFs_test  = list(CM_TFs, ENDO1_TFs, FIB1_TFs, FIB2_TFs, MACROPH_TFs)

dorothea_mat = as.matrix(visium_slide[["dorothea"]]@data)


healthy_results = pmap(all_res[,c("cell_id","TF_act_matB","TFs_test")], function(cell_id,TF_act_matB,TFs_test){
  
  path = paste0("./visium_results_manuscript/atac_validation/healthy/", cell_id,"_modulescore.pdf")
  
  modules_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                            TF_act_matB = as.data.frame(TF_act_matB),
                            plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  
  path = paste0("./visium_results_manuscript/atac_validation/healthy/", cell_id,"_dorothea.pdf")
  

  if(sum(TFs_test %in% rownames(dorothea_mat))>0){
    dorothea_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                               TF_act_matB = as.data.frame(dorothea_mat),
                               plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  }else{
    dorothea_res = NULL
  }

  return(enframe(list("module_score" = modules_res,
                      "dorothea_score" = dorothea_res)) %>% unnest())
})

names(healthy_results) = all_res$cell_id
write.table(file = "./visium_results_manuscript/atac_validation/healthy/wilcoxon_results.txt",
            enframe(healthy_results) %>% unnest(),
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")


# Specific analysis throughout the paper
# Infarct - IZ
atac = "CK174"
slide = "157775"

# Read visium
slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                     slide,slide)

visium_slide = readRDS(slide_file)
DefaultAssay(visium_slide) = "SCT"

# Getting the regulatory networks
atac_path = sprintf("./MotifToGene/%s",atac)
all_res = readRDS(file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac))

TFs = c("NRF1","MEF2C")

# Generating the MS matrices per cell-type

all_res = all_res %>%
  dplyr::filter(cell_id %in% c("Cardiomyocyes","Fibroblass_0")) %>%
  group_by(cell_id) %>%
  dplyr::mutate(TF_act_matB = map(MS_regulon, function(x){
    
    TF_act = getTF_matrix_MS(visium_slide,MS_regulon = x[TFs])
    
    return(TF_act)
  }))

# Testing them

all_res$TFs_test  = list(TFs, TFs)

dorothea_mat = as.matrix(visium_slide[["dorothea"]]@data)
hint_mat = as.matrix(visium_slide[["HINT_TF_Activity"]]@data)

ischemic_results = pmap(all_res[,c("cell_id","TF_act_matB","TFs_test")], function(cell_id,TF_act_matB,TFs_test){
  
  path = paste0("./visium_results_manuscript/atac_validation/ischemic/", cell_id,"_modulescore.pdf")
  
  modules_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                            TF_act_matB = as.data.frame(TF_act_matB),
                            plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  
  path = paste0("./visium_results_manuscript/atac_validation/ischemic/", cell_id,"_dorothea.pdf")
  
  if(sum(TFs_test %in% rownames(dorothea_mat))>0){
    dorothea_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                               TF_act_matB = as.data.frame(dorothea_mat),
                               plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  }else{
    dorothea_res = NULL
  }
  
  pdf(file = paste0("./visium_results_manuscript/atac_validation/ischemic/", cell_id,"_spatialplots.pdf"))
  
  DefaultAssay(visium_slide) = "HINT_TF_Activity"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  visium_slide[['TF_act_matB_MS']] = CreateAssayObject(data = TF_act_matB)
  DefaultAssay(visium_slide) = "TF_act_matB_MS"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  DefaultAssay(visium_slide) = "dorothea"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  dev.off()
  
  return(enframe(list("module_score" = modules_res,
                      "dorothea_score" = dorothea_res)) %>% unnest())
})

names(ischemic_results) = all_res$cell_id
write.table(file = "./visium_results_manuscript/atac_validation/ischemic/wilcoxon_results.txt",
            enframe(ischemic_results) %>% unnest(),
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# Specific analysis throughout the paper
# Border Zone
atac = "CK168"
slide = "157781"

# Read visium
slide_file = sprintf("./visium_results_manuscript/processed_objects/%s.rds",
                     slide,slide)
visium_slide = readRDS(slide_file)
DefaultAssay(visium_slide) = "SCT"

# Getting the regulatory networks
atac_path = sprintf("./MotifToGene/%s",atac)
all_res = readRDS(file = sprintf(paste0(atac_path,"/%s_TFvalidation.rds"),atac)) %>%
  mutate(MS_regulon = map(MS_regulon, function(x){
    names(x) = toupper(names(x))
    return(x)
  }))

CM_TFs = c("GATA4", "MEF2C","MYOD1")
CM1_TFs = c("NFE2L1")
CM2_TFs = c("SMAD3","SMAD4","SMAD5", "JUN", "FOS")
FB_TFs = c("RUNX1","RUNX2")

TFs = c(CM_TFs, CM1_TFs, CM2_TFs, FB_TFs)

all_res = all_res %>%
  dplyr::filter(cell_id %in% c("Cardiomyocyes_1",
                               "Cardiomyocyes_2",
                               "Fibroblass_1")) %>%
  group_by(cell_id) %>%
  dplyr::mutate(TF_act_matB = map(MS_regulon, function(x){
    
    TF_act = getTF_matrix_MS(visium_slide,MS_regulon = x[TFs])
    
    return(TF_act)
  }))

all_res$TFs_test  = list(TFs, TFs, TFs)

dorothea_mat = as.matrix(visium_slide[["dorothea"]]@data)
hint_mat = as.matrix(visium_slide[["HINT_TF_Activity"]]@data)

bz_results = pmap(all_res[,c("cell_id","TF_act_matB","TFs_test")], function(cell_id,TF_act_matB,TFs_test){
  
  path = paste0("./visium_results_manuscript/atac_validation/borderzone/", cell_id,"_modulescore.pdf")
  
  modules_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                            TF_act_matB = as.data.frame(TF_act_matB),
                            plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  
  path = paste0("./visium_results_manuscript/atac_validation/borderzone/", cell_id,"_dorothea.pdf")
  
  
  if(sum(TFs_test %in% rownames(dorothea_mat))>0){
    dorothea_res = (compareTFs(TF_act_matA = as.data.frame(hint_mat),
                               TF_act_matB = as.data.frame(dorothea_mat),
                               plotting = T,features = TFs_test, quantileA = 0.90, path = path))
  }else{
    dorothea_res = NULL
  }
  
  pdf(file = paste0("./visium_results_manuscript/atac_validation/borderzone/", cell_id,"_spatialplots.pdf"))
  
  DefaultAssay(visium_slide) = "HINT_TF_Activity"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  visium_slide[['TF_act_matB_MS']] = CreateAssayObject(data = TF_act_matB)
  DefaultAssay(visium_slide) = "TF_act_matB_MS"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  DefaultAssay(visium_slide) = "dorothea"
  plot(SpatialFeaturePlot(visium_slide, features = TFs_test))
  
  dev.off()
  
  return(enframe(list("module_score" = modules_res,
                      "dorothea_score" = dorothea_res)) %>% unnest())
})

names(bz_results) = all_res$cell_id
write.table(file = "./visium_results_manuscript/atac_validation/borderzone/wilcoxon_results.txt",
            enframe(bz_results) %>% unnest(),
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
