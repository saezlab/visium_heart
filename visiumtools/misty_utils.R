# Copyright (c) [2020] [Ricardo O. Ramirez Flores, Jovan Tanevski]
# roramirezf@uni-heidelberg.de

#' Utilities to plot and analyse MISTy outputs
#' 
#' Most of the code is from Jovan adapted for general usage
#' 
#' 

#' Function to get optimal results
#' Here, first we identify which l parameter is the best for each target
#' Then we create an optimal directory and copy specific information
#' 
#' @param out_dir_name = out_alias used in MISTy runs
#' @param ls = l parameters used in the MISTy runs
#' @return a folder out_dir_name_optim with optimal info
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
  
  ## to do: save optimal results
  
  write.table(optimal.l,
              file = paste0(out_dir_name,"_optim/optim_l.txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  
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

#'Function to get processed MISTys final output
#'@param results_folder: result folder from a MISTy run or an optimization search
#'@param p.cutoff: p-value cutoff to evaluate coefficients
#'
#'@return a list with improvement, contribution, coefficients and importance info
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
    purrr::reduce(`+`)) / length(images)) %>% 
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
        dplyr::select((order(colnames(result))))
      #dplyr::select(noquote(order(colnames(result))))
    })
  })
  
  aggregated = maps %>% purrr::reduce(function(acc, l) {
    map2(acc, l, ~ (((.x %>% dplyr::select(-Predictor)) + (.y %>% dplyr::select(-Predictor))) %>%
                      mutate(Predictor = .x %>% pull(Predictor))))
  })
  
  performance_df = read_delim(paste0(results_folder,"/performance.txt"), delim = " ")
  
  result_list = list("impr" = impr,
                     "targets" = targets,
                     "coefs" = coefs,
                     "importance" = aggregated,
                     "performance" = performance_df)
  
  return(result_list)
  
}


#'Plot MISTy performance
#'@param MISTy_out: result list from MISTy_aggregator
#'@param pdf_out: pdf file where to return plots
#'@param predicted_features: show results of selected features, if NULL all
#'
#'@return a pdf with performance plots

plot_misty_performance = function(MISTy_out,
                                  pdf_out,
                                  predicted_features = NULL){
  #Generating figures
  #Overall performance
  performance_df = MISTy_out$performance %>% tidyr::pivot_longer(cols = -target) %>%
    dplyr::filter(grepl("R2",name))
  
  R2_impr = MISTy_out$impr %>% dplyr::select(contains("R2")) %>%
    mutate(target = targets) %>%
    pivot_longer(cols = -target, 
                 names_to = "name", 
                 values_to = "value") %>% 
    arrange(desc(value))
  
  #Look at importances of views
  contribution_df = MISTy_out$coefs
  
  if(length(predicted_features)>0){
    
    performance_df = performance_df %>%
      dplyr::filter(target %in% predicted_features)
    
    R2_impr = R2_impr %>%
      dplyr::filter(target %in% predicted_features)
    
    contribution_df = contribution_df %>%
      dplyr::filter(target %in% predicted_features)
    
  }
  
  performance_plt = ggplot(performance_df,
                             aes(fill = name, y = target, x = value)) +
    geom_bar(stat = "identity",position="dodge") + 
    theme_minimal()
  
  impr_plot = ggplot(R2_impr) +
    geom_point(aes(x = target, y = value * 100)) +
    theme_classic() +
    ylab("Change in variance explained") +
    xlab("Target") +
    theme(axis.title = element_text(size=11),
          axis.text = element_text(size=10),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  coefs_plot = ggplot(contribution_df) + 
    geom_col(aes(x=target, 
                 y=value, group=view, fill=view)) +
    xlab("Target") +
    ylab("Contribution") +
    theme_classic() +
    theme(axis.text = element_text(size=13),
          axis.title = element_text(size=13),
          legend.text = element_text(size=11)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                     vjust = 0.5)) +
    scale_fill_brewer(palette="Dark2",
                      labels = c("Intrinsic","Para pathway")) +
    labs(fill = "View")
  
  pdf(file = pdf_out, width = 12, height = 9,onefile = TRUE)
  plot(performance_plt)
  plot(impr_plot)
  plot(coefs_plot)
  dev.off()
}

#'Plot view importances: To do, expand to multiple views with a walk through views
#'@param MISTy_out: result list from MISTy_aggregator
#'@param pdf_out: pdf file where to return plots
#'@param importance_cut: to filter predictors of no use
#'@return a pdf document with importance plots
plot_misty_importance = function(MISTy_out,
                                 pdf_out,
                                 predicted_features = NULL,
                                 predictors_features = NULL,
                                 importance_cut){
  
  importance_intra = tidyr::gather(MISTy_out$importance[[1]], 
                                   "Predicted",
                                   "Importance", -Predictor)
  
  importance_para = tidyr::gather(MISTy_out$importance[[2]], 
                                  "Predicted",
                                  "Importance", -Predictor)
  
  if(length(predicted_features)>0){
    
    importance_intra = importance_intra %>% 
      dplyr::filter(Predicted %in% predicted_features)
    
    importance_para =  importance_para %>% 
      dplyr::filter(Predicted %in% predicted_features)
  }
  
  if(length(predictors_features)>0){
    
    importance_intra = importance_intra %>% 
      dplyr::filter(Predictor %in% predictors_features)
    
    importance_para =  importance_para %>% 
      dplyr::filter(Predictor %in% predictors_features)
  }
  
  #INTRA

  intra_summ = importance_intra %>% 
    dplyr::mutate(importance_bool = Importance >= importance_cut) %>%
    group_by(Predictor) %>% 
    summarize(predictor_summ = sum(importance_bool,na.rm = T)) %>%
    dplyr::filter(predictor_summ >= 1)
  
  intra_predictors = intra_summ %>% select(Predictor) %>% pull()
  
  importance_intra = importance_intra %>% 
    dplyr::filter(Predictor %in% intra_predictors) %>%
    dplyr::arrange(Predicted, -Importance)
  
  intra_plot = importance_intra %>% 
    ggplot(aes(x = factor(Predictor,
                          levels = unique(importance_intra$Predictor)), 
               y = factor(Predicted,
                          levels = unique(importance_intra$Predicted)) ,
               fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text = element_text(size=10),
          legend.key.size = unit(.6, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "white", 
                         mid = "white", 
                         high = scales::muted("blue"),
                         midpoint = 0.9) +
    xlab("Intrinsic features")
  
  #PARA
  para_summ = importance_para %>% 
    dplyr::mutate(importance_bool = Importance >= importance_cut) %>%
    group_by(Predictor) %>% 
    summarize(predictor_summ = sum(importance_bool,na.rm = T)) %>%
    dplyr::filter(predictor_summ >= 1)
  
  para_predictors = para_summ %>% select(Predictor) %>% pull()
  
  importance_para = importance_para %>%
    dplyr::arrange(Predicted, -Importance) %>%
    dplyr::filter(Predictor %in% para_predictors)
  
  para_plot = importance_para %>% 
    ggplot(
      aes(x = factor(Predictor,
                     levels = unique(importance_para$Predictor)), 
          y = factor(Predicted,
                     levels = unique(importance_para$Predicted)),
          fill = Importance)) + geom_tile() + 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size=11),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text = element_text(size=10),
          legend.key.size = unit(.6, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.position = "bottom") +
    scale_fill_gradient2(low = "white", 
                         mid = "white", 
                         high = scales::muted("blue"),
                         midpoint = 0.9) +
    xlab("Para features")
  
  pdf(file = pdf_out, width = 12, height = 13,onefile = TRUE)
  plot(intra_plot)
  plot(para_plot)
  dev.off()
  
  importance_out = list("intra"= importance_intra,
                        "para" = importance_para)
  
}

#'Gets significant predictions from MISTy
#'@param results_folder: folder path containing MISTy runs
#'@param p_value_trsh: to filter significant changes in importance
#'@param R2_trsh: to filter fits greater than the stated R
#'@return 2 pdf documents with performance and importance plots
get_misty_bic = function(results_folder,
                   p_value_trsh = 0.15,
                   R2_trsh = 0){
  
  MISTyout = MISTy_aggregator(results_folder = results_folder,
                              p.cutoff = 0.05)
  
  bic = MISTyout$performance %>% 
    dplyr::arrange(p.R2) %>% 
    dplyr::filter(p.R2 <= p_value_trsh,
                  intra.R2 >= R2_trsh) %>% 
    select(target) %>% pull()
  
  return(bic)
}



#'Gets significant predictions from MISTy
#'@param MISTy_out_folders: vector with folder paths containing MISTy runs
#'@param p_value_trsh: to filter significant changes in importance
#'@param R2_trsh: to filter fits greater than the stated R
#'@param performance_alias: output alias
#'@param importance_alias: output alias
#'@return 2 pdf documents with performance and importance plots
plot_misty_bic = function(MISTy_out_folders,
                          p_value_trsh = 0.15,
                          R2_trsh = 0,
                          performance_alias = "_prfmnce.pdf",
                          importance_alias = "_imp.pdf",
                          importance_cut = 1){
  
  for(f in MISTy_out_folders){
    
    MISTyout = MISTy_aggregator(results_folder = f,
                                p.cutoff = 0.05)
    
    bic = MISTyout$performance %>% 
      dplyr::arrange(p.R2) %>% 
      dplyr::filter(p.R2 <= p_value_trsh,
                    intra.R2 >= R2_trsh) %>% 
      select(target) %>% pull()
    
    if(length(bic) > 0){
    
    plot_misty_performance(MISTy_out = MISTyout,
                           predicted_features = bic,
                           pdf_out = paste0(f,
                                            performance_alias))
    
    plot_misty_importance(MISTy_out = MISTyout,
                          pdf_out = paste0(f,
                                           importance_alias),
                          predicted_features = bic,
                          predictors_features = NULL,
                          importance_cut = importance_cut)
    
    }
    
  }
}



#'Gets significant predictions from MISTy
#'@param MISTy_out: R object coming from MISTy_aggregator
#'@param feature_ann: to filter significant changes in importance
#'@return MISTy_out object with modified importances
group_MISTy_out = function(MISTy_out,
                          feature_ann){
  
  #MISTy_out
  MISTy_out_group = MISTy_out
  
  importance_intra = tidyr::gather(MISTy_out$importance[[1]], 
                                   "Predicted",
                                   "Importance", -Predictor) %>%
    dplyr::filter(Predicted %in% feature_ann$feature) %>%
    left_join(feature_ann, by = c("Predicted" = "feature")) %>%
    group_by(Predictor,ann) %>% 
    summarise(Importance = mean(Importance,na.rm = T)) %>%
    spread(key = ann,value = Importance)
    
  
  importance_para = tidyr::gather(MISTy_out$importance[[2]], 
                                  "Predicted",
                                  "Importance", -Predictor) %>%
    dplyr::filter(Predicted %in% feature_ann$feature) %>%
    left_join(feature_ann, by = c("Predicted" = "feature")) %>%
    group_by(Predictor,ann) %>% 
    summarise(Importance = mean(Importance,na.rm = T)) %>%
    spread(key = ann,value = Importance)
  
  MISTy_out_group$importance[[1]] = importance_intra
  MISTy_out_group$importance[[2]] = importance_para
  
  return(MISTy_out_group)
}

#'Gets para expression from MISTy
#'@param visium_slide: Seurat object with assays to be transformed
#'@param para_assay: slot in visium slide to be transformed
#'@param para_features: para features to transform
#'@param l: radius cost parameter
#'@return a matrix ready to be included as an assay in the visium object
get_para_matrix = function(visium_slide, para_assay, para_features, l = 10){
  clear_cache()
  # Getting data ready to create views
  geometry = visium_slide@images$slice1@coordinates[,c(2,3)]
  
  # Para data
  para_df = as.matrix(visium_slide@assays[[para_assay]]@data)
  para_df = para_df %>%
    t %>% data.frame(check.names = F)
  para_df = para_df[rownames(geometry),]
  
  # Defining useful data para
  para_df = para_df[rownames(geometry),para_features]
  
  colnames(para_df) = gsub("-","_", colnames(para_df))
  
  views_para = create_initial_view(para_df, 
                                   unique.id = paste0("para_",l^2)) %>% 
    add_paraview(geometry, l^2)
  
  # Fetching actual para info to be used
  
  # Spot specific view comes from the view above
  data_red = views_para[[3]]$data
  rownames(data_red) = rownames(views_para$intracellular$data) #we named rows just for easy access
  
  para_mat = (t(as.matrix(data_red)))[,colnames(visium_slide)]
  
  return(para_mat)
}







