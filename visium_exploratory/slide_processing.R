# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Processing of breast cancer visium slides
#' It normalizes the data using SCT transform,
#' calculates TF and PROGENy activities and
#' defines a matrix of expressed ligands
#' 

library(Seurat)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(viper)
library(progeny)
library(OmnipathR)
source("MISTy-master/multiview.R")


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

## Function to process SpaceRanger folders
## Input: Folder path
## Output: Seurat object ready for MISTy
process_visium = function(dir_path){
  ## <<Gene_expression>>
  dat = Load10X_Spatial(data.dir = dir_path)
  dat = SCTransform(dat, assay = "Spatial", verbose = FALSE)
  
  ## <<DOROTHEA>>
  dorothea_regulon_human = read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")
  # We obtain the regulons based on interactions with confidence level A, B and C
  regulon = dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%
    df2regulon()
  
  tf_act_mat = viper(eset = as.matrix(dat[["SCT"]]@data), 
                     regulon = regulon, nes = TRUE, 
                     method = "scale", minsize = 4, 
                     eset.filter = FALSE,
                     verbose = FALSE)
  
  ## Repeated regulon RFXAP
  tf_act_mat = tf_act_mat[!duplicated(tf_act_mat),]
  dat[['dorothea']] = CreateAssayObject(counts = tf_act_mat)
  
  
  ## <<PROGENY>>
  progeny_scores = progeny::progeny(expr = as.matrix(dat[["SCT"]]@data),
                                    scale=TRUE, 
                                    organism="Human", top=1000, perm=1)
  
  dat[['progeny']] = CreateAssayObject(counts = t(progeny_scores))

  ## <<Ligands>>
  
  opath_ligands = readRDS(file = "./markers/opath_ligands.rds")
  ligands = opath_ligands[opath_ligands$GeneSymbol %in% 
                          rownames(dat@assays$SCT@data),"GeneSymbol"] 
  
  ligands_gex = as.matrix(dat@assays$SCT@data[ligands,])
  ligands_expr = ligands_gex[rowSums(ligands_gex > 0)/ncol(ligands_gex) > .3,]
  
  dat[['ligands']] = CreateAssayObject(counts = ligands_expr)
  
  ## <<Heart markers scores>>
  markers_mat = readRDS("./markers/markers_matrix.rds")
  exprdata = as.matrix(dat[["SCT"]]@data)
  genes_in_markers = rownames(exprdata) %in% rownames(markers_mat)
  exprdata = exprdata[genes_in_markers,]
  markers_mat = markers_mat[rownames(exprdata),]
  
  ct_scores = scale(t(t(markers_mat) %*% exprdata))
  
  dat[['ctscores']] = CreateAssayObject(counts = t(ct_scores))
  
  return(dat)
}

## Function to run MISTy using ligands as paracrine view and 
## PROGENy scores as target features

run_ligand_MISTy = function(seurat_visium_obj,
                            ligs_features,
                            l = 20,
                            out_dir_name = "default"){
  
  plan(multiprocess, workers = 6)
  
  # Getting data ready to create views
  
  geometry = seurat_visium_obj@images$slice1@coordinates
  
  pth = as.matrix(seurat_visium_obj@assays$progeny@data) %>% 
    t %>% data.frame 
  
  pth = pth[rownames(geometry),]
  
  ligs = as.matrix(seurat_visium_obj@assays$ligands@data) %>% 
    t %>% data.frame(check.names = F)
  
  ligs = ligs[rownames(geometry),ligs_features]
  
  # Creating views
  
  views_pth = create_initial_view(pth, unique.id = paste0(out_dir_name,"path","_",l^2)) %>% 
    add_paracrine_view(geometry[ ,2:3], l^2)
  
  views_ligs = create_initial_view(ligs, unique.id = paste0(out_dir_name,"ligs","_",l^2)) %>% 
    add_paracrine_view(geometry[ ,2:3], l^2)
  
  # Here we define the paraview of ligands
  views = views_pth %>% 
    add_views(create_view(paste0(out_dir_name,"ligs_para_",l^2),
                          views_ligs[[3]]$data))
  
  MISTy_run = estimate_importances(views, 
                                   paste0(out_dir_name,"_",l^2))
}

## Function to get optimal results
## Here, first we identify which l parameter is the best for each target
## Then we create an optimal directory and copy specific information
## Inputs:
## out_dir_name: directory name used in run_ligand_MISTy
## ls: l's tried in run_ligand_MISTy

get_optimal = function(out_dir_name,
                       ls){
  
  system(paste0("mkdir ", out_dir_name,"_optim"))
  
  l = ls
  
  perf = (l^2) %>% map_dfr(function(p){
    performance <- read_delim(paste0(out_dir_name, "_",p, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    #selection criterion maximum R2 performance improvement per marker
    performance %>% arrange(target) %>%
      mutate(impr = multi.R2 - intra.R2) %>% 
      dplyr::select(target,impr) %>% 
      column_to_rownames("target") %>% t %>% as.data.frame
  })
  
  # For each target get the maximum of improvement
  # distinct ls
  optimal.l = colnames(perf) %>% enframe(name = NULL, value = "target") %>% 
    mutate(l = apply(perf, 2, which.max) %>% (l^2)[.])
  
  # Copy the relevant documents to the new directory, deleting the l parameter info
  optimal.l %>% pull(target) %>% walk2(optimal.l %>% pull(l), function(.x, .y){
    files <- list.files(paste0(out_dir_name,"_", .y, "/"), 
                        paste0("importances_", .x, '*'), 
                        full.names = TRUE)
    
    files %>% walk(~file.copy(., paste0(out_dir_name,"_optim/", 
                                        str_replace(last(str_split(., "/")[[1]]), 
                                                    "a([\\._][0-9]+)+", "a"))))
  })
  
  #very suboptimal
  optimal.l %>% pull(target) %>% walk2(optimal.l %>% pull(l), function(.x, .y){
    performance <- read_delim(paste0(out_dir_name,"_", .y, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    coeff <- read_delim(paste0(out_dir_name,"_", .y, "/coefficients.txt"),
                        delim = " ", col_types = cols()
    ) %>% distinct()
    
    
    
    if(!file.exists(paste0(out_dir_name,"_optim/performance.txt"))){
      write(colnames(performance) %>% paste0(collapse=" "), 
            paste0(out_dir_name,"_optim/performance.txt"))
    }
    
    if(!file.exists(paste0(out_dir_name,"_optim/coefficients.txt"))){
      write(str_remove_all(colnames(coeff) %>% paste0(collapse=" "), "([\\._][0-9]+)+"), 
            paste0(out_dir_name,"_optim/coefficients.txt"))
    }
    
    write(performance %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/performance.txt"), append = T)
    
    write(coeff %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/coefficients.txt"), append = T)
    
  })
  
}

## Function to get MISTys final output
## Inputs:
## results_folder: a vector with directory names of optimal results
## p.cutoff: p-value cutoff
## Outputs:
## a list with improvement, contribution, coefficients and importance info

MISTy_aggregator = function(results_folder,
                            p.cutoff = 0.05){
  
  images = results_folder
  
  # Improvement 
  impr = images %>% map_dfc(function(image) {
    performance <- read_delim(paste0(image, .Platform$file.sep,
                                     "performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    targets <<- unique(performance$target)
    performance %>%
      arrange(target) %>%
      transmute(RMSE = (intra.RMSE - multi.RMSE) / intra.RMSE, R2
                = (multi.R2 - intra.R2))
  })
  
  impr_plot = ggplot(impr %>% dplyr::select(contains("RMSE")) %>%
           mutate(target = targets) %>%
           tidyr::pivot_longer(cols = -target, names_to = "name", values_to = "value")) +
    geom_boxplot(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Improvement (%)") +
    xlab("Target") +
    ylim(c(-5, 17)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  avg <- ((images %>% map(function(image) {
    coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
                               delim = " ", col_types = cols()
    )
    
    targets <<- coefficients %>%
      pull(target) %>%
      sort
    
    coefficients %>%
      distinct() %>%
      arrange(target) %>%
      dplyr::select(-target, -contains("intercept")) %>%
      mutate_at(vars(starts_with("p.")), ~ as.numeric(. <= p.cutoff)) %>%
      mutate_at(vars(-starts_with("p.")), abs)
  }) %>% 
    reduce(`+`)) / length(images)) %>% 
    mutate(target = targets)
  
  ctotals <- avg %>% 
    dplyr::select(-starts_with("p."), -"target") %>%
    rowSums
  
  coefs <- avg %>% 
    dplyr::select(-starts_with("p.")) %>%
    mutate_if(is.numeric, ~./ctotals) %>%
    tidyr::pivot_longer(-target, names_to = "view")
  
  maps = images %>% map(function(image) {
    coefficients <- read_delim(paste0(image, .Platform$file.sep, "coefficients.txt"),
                               delim = " ", col_types = cols()
    ) %>% distinct()
    
    targets <- unique(coefficients$target)
    views <<- (coefficients %>% dplyr::select(-target, -starts_with("p."), -intercept) %>% colnames())
    
    # one heatmap per view
    maps <- views %>% map(function(view) {
      all.importances <- targets %>% map(~ read_csv(paste0(
        image, .Platform$file.sep, "importances_",
        .x, "_", view, ".txt"
      ),
      col_types = cols()
      ) %>%
        distinct() %>%
        filter(!grepl("_2$", target)))
      
      features <- unique(all.importances %>% map(~ .x$target) %>% unlist())
      
      pview <- paste0("p.", view)
      ps <- coefficients %>%
        dplyr::select(target, !!pview) %>%
        mutate(!!pview := (1 - !!sym(pview)))
      
      
      # importances are standardized for each target an multiplied by 1-pval(view)
      result <- all.importances %>%
        imap_dfc(~
                   tibble(target = features, zero.imp = 0) %>%
                   left_join(.x, by = "target") %>%
                   transmute(feature = target, importance = (zero.imp + scale(imp)[, 1]) *
                               (ps %>% filter(target == targets[.y]) %>% pull(pview))) %>%
                   dplyr::select(importance)) %>%
        `colnames<-`(targets) %>%
        mutate(Predictor = features)
      
      # in order for aggregation
      result %>%
        arrange(Predictor) %>%
        dplyr::select(noquote(order(colnames(result))))
    })
  })
  
  aggregated = maps %>% reduce(function(acc, l) {
    map2(acc, l, ~ (((.x %>% dplyr::select(-Predictor)) + (.y %>% dplyr::select(-Predictor))) %>%
                      mutate(Predictor = .x %>% pull(Predictor))))
  })
  
  result_list = list("impr" = impr,
                     "targets" = targets,
                     "coefs" = coefs,
                     "importance" = aggregated)
  
  return(result_list)
  
}



