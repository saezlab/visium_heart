# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Wrappers of niche characterizers

library(compositions)
library(uwot)
library(factoextra)
library(ggpubr)

annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

patient_time <- read_csv("./markers/visium_timeinfo.csv") %>%
  dplyr::select(sample_id, `days after infarction`) %>%
  dplyr::rename("time" = `days after infarction`) %>%
  dplyr::mutate(time = ifelse(time == "control", 0, time)) %>%
  dplyr::mutate(time = as.numeric(time))

patient_info <- read_csv("./markers/visium_patient_anns_revisions.csv") %>%
  left_join(annotation_names) %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))) %>%
  left_join(patient_time)

# 1. Get the proportions of the niches

get_niche_props <- function(niche_info) {
  
  niche_info %>%
    group_by(patient_region_id, mol_niche) %>%
    summarize(n_spots = n()) %>%
    mutate(n_spots_pat = sum(n_spots)) %>%
    mutate(niche_prop = n_spots/n_spots_pat) %>%
    ungroup()
  
}

# 2. Generate plots based on proportions -----------------------------------------------------------

plot_dots_niche <- function(niche_props) {
  
  dot_plt <- ggplot(niche_props, aes(x = mol_niche, y = patient_region_id, size = niche_prop)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.y = element_text(size=12)) +
    xlab("") +
    ylab("")
  
}

# Do patient group comparisons

# Filter compositions of niches per sample 
# Filter slides whenever a niche represents more than 80% of spots
# Complete information with 0s

filter_compositions <- function(niche_props, by = "patient") {
  
  if(by == "patient") {
    
    high_qc_pats <- niche_props %>%
      group_by(patient_region_id) %>%
      summarize(max_comp =  max(niche_prop)) %>%
      dplyr::filter(max_comp < 0.8) %>%
      pull(patient_region_id)
    
    filtered_niche_props <- niche_props %>%
      dplyr::filter(patient_region_id %in% high_qc_pats) %>%
      dplyr::select(patient_region_id, mol_niche, niche_prop) %>%
      tidyr::complete(patient_region_id, mol_niche, fill = list("niche_prop" = 0)) %>%
      left_join(dplyr::select(patient_info, patient_region_id, patient_group, major_labl, time) %>% unique())
    
  } else if(by == "niche") {
    
    high_qc_niches <- niche_props %>%
      dplyr::filter(niche_prop > 0.01) %>%
      dplyr::group_by(mol_niche) %>%
      summarise(n_pats = n()) %>%
      dplyr::filter(n_pats >= 5) %>%
      pull(mol_niche)
    
    filtered_niche_props <- niche_props %>%
      dplyr::filter(mol_niche %in% high_qc_niches) %>%
      dplyr::select(patient_region_id, mol_niche, niche_prop) %>%
      tidyr::complete(patient_region_id, mol_niche, fill = list("niche_prop" = 0)) %>%
      left_join(dplyr::select(patient_info, patient_region_id, patient_group, major_labl, time) %>% unique())
    
    
  }
  
 
  
}

# Generate box_plots

plot_box_niches <- function(filtered_niche_props) {
  
  box_plt <- ggplot(filtered_niche_props, aes(x = patient_group, y = niche_prop, color = patient_group)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
    facet_wrap(.~ mol_niche, ncol = 4, scales = "free") +
    ylab("niche proportion")
 
   
}

plot_box_niches_pw <- function(filtered_niche_props, pw_niche_prop_comps, area = F) {
  
  props_info <- filtered_niche_props %>%
    group_by(mol_niche) %>%
    nest() %>%
    dplyr::rename("props_data" = data)
    
  pw_info <- pw_niche_prop_comps %>%
      group_by(mol_niche) %>%
      dplyr::select(group1, group2, p , p.adj) %>%
      nest() %>%
      dplyr::rename("pw_data" = data)
  
  boxplot_df <- left_join(props_info, pw_info) %>%
    dplyr::filter(map_lgl(pw_data, ~ is.null(.x)) != TRUE) %>%
    mutate(bxplt = map2(props_data, pw_data, function(props, pw, area_flag = area) {
      
      # Here you filter areas that were tested
      if(area_flag) {
        groups_list <- unique(c(pw$group1, pw$group2))
        props <- props %>%
          dplyr::filter(major_labl %in% groups_list)
      }
      
      max_val <- max(props$niche_prop) + 0.05
      
      if(area_flag) {
        
        box_plt <- ggplot(props, aes(x = major_labl, y = niche_prop, color = major_labl))
        
      } else {
        
        box_plt <- ggplot(props, aes(x = patient_group, y = niche_prop, color = patient_group))
        
      }
  
      box_plt <- box_plt +
        geom_boxplot() +
        geom_point() +
        ggpubr::stat_pvalue_manual(pw, label = "p.adj", 
                                   y.position = max_val, 
                                   step.increase = 0.1,
                                   tip.length = 0.01,size = 3) +
        theme_classic() +
        theme(axis.text = element_text(size = 11),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.5, 
                                         size = 11),
              legend.position = "none") +
        ylab("niche proportion") 
      
      
    }))
  
  return(boxplot_df)
}

plot_box_niches_pw_area <- function(filtered_niche_props, pw_niche_prop_comps) {
  
  props_info <- filtered_niche_props %>%
    group_by(mol_niche) %>%
    nest() %>%
    dplyr::rename("props_data" = data)
  
  pw_info <- pw_niche_prop_comps %>%
    group_by(mol_niche) %>%
    dplyr::select(group1, group2, p , p.adj) %>%
    nest() %>%
    dplyr::rename("pw_data" = data)
  
  boxplot_df <- left_join(props_info, pw_info) %>%
    mutate(bxplt = map2(props_data, pw_data, function(props, pw) {
      
      max_val <- max(props$niche_prop) + 0.05
      
      box_plt <- ggplot(props, aes(x = patient_group, y = niche_prop, color = patient_group)) +
        geom_boxplot() +
        geom_point() +
        ggpubr::stat_pvalue_manual(pw, label = "p.adj", 
                                   y.position = max_val, 
                                   step.increase = 0.1,
                                   tip.length = 0.01,size = 3) +
        theme_classic() +
        theme(axis.text = element_text(size = 11),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.5, 
                                         size = 11),
              legend.position = "none") +
        ylab("niche proportion") 
      
      
    }))
  
  return(boxplot_df)
}


plot_box_niches_area <- function(filtered_niche_props, myogenic_enriched = F) {
  
  if(myogenic_enriched) {
    filtered_niche_props <- filtered_niche_props %>%
      dplyr::filter(patient_group == "myogenic-enriched",
                    ! major_labl %in% c("FZ", "IZ"))
  } 
  
  box_plt <- ggplot(filtered_niche_props, aes(x = major_labl, y = niche_prop, color = major_labl)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
    facet_wrap(.~ mol_niche, ncol = 4, scales = "free") +
    ylab("niche proportion")
  
  
}

# Calculate Kruskall-Wallis

kw_niche_prop_test <- function(filtered_niche_props) {
  
  niche_test_kw <- filtered_niche_props %>%
    dplyr::select(mol_niche, niche_prop, patient_group) %>%
    group_by(mol_niche) %>%
    nest() %>%
    mutate(kw_res = map(data, function(dat) {
      
      kruskal.test(niche_prop ~ patient_group, 
                   data = dat) %>%
        broom::tidy()
      
      
    })) %>%
    dplyr::select(mol_niche, kw_res) %>%
    unnest() %>%
    ungroup() %>%
    mutate(corr_pval = p.adjust(p.value))
  
}

kw_niche_prop_test_area <- function(filtered_niche_props, myogenic_enriched = F) {
  
  if(myogenic_enriched) {
    filtered_niche_props <- filtered_niche_props %>%
      dplyr::filter(patient_group == "myogenic-enriched",
                    ! major_labl %in% c("FZ", "IZ"))
  } 
  
  niche_test_kw <- filtered_niche_props %>%
      dplyr::select(mol_niche, niche_prop, major_labl) %>%
      group_by(mol_niche) %>%
      nest() %>%
      mutate(kw_res = map(data, function(dat) {
        
        kruskal.test(niche_prop ~ major_labl, 
                     data = dat) %>%
          broom::tidy()
        
        
      })) %>%
      dplyr::select(mol_niche, kw_res) %>%
      unnest() %>%
      ungroup() %>%
      mutate(corr_pval = p.adjust(p.value))
  
}

# Pairwise comparisons

pw_niche_prop_test <- function(filtered_niche_props,
                               area = F,
                               myogenic_enriched = F) {
  
  if(area) {
    
    if(myogenic_enriched) {
      filtered_niche_props <- filtered_niche_props %>%
        dplyr::filter(patient_group == "myogenic-enriched",
                      ! major_labl %in% c("FZ", "IZ"))
    } 
    
    return(filtered_niche_props %>%
      group_by(mol_niche) %>%
      nest() %>%
      dplyr::mutate(pw_res = map(data, function(dat) {
        
        compare_means(niche_prop ~ major_labl,  
                      data = dat,
                      method = "wilcox.test", 
                      alternative = "two.sided")
        
        
      })) %>%
      dplyr::select(pw_res) %>%
      unnest())
    
  } else {
    
    return(filtered_niche_props %>%
      group_by(mol_niche) %>%
      nest() %>%
      dplyr::mutate(pw_res = map(data, function(dat) {
        
        compare_means(niche_prop ~ patient_group,  
                      data = dat,
                      method = "wilcox.test", 
                      alternative = "two.sided")
        
        
      })) %>%
      dplyr::select(pw_res) %>%
      unnest())
    
    
  }
}


# 3. ILR transformations, clustering, UMAPs and PCAs

nichedf_tomatrix <- function(filtered_niche_props) {
  
  complete_niche_info_mat <- filtered_niche_props %>%
    dplyr::select(-c("patient_group", "time", "major_labl")) %>%
    pivot_wider(names_from = mol_niche, values_from = niche_prop) %>%
    as.data.frame() %>%
    column_to_rownames("patient_region_id") %>%
    as.matrix()
  
}

# columns niches
# rows patients

ILR_transform <- function(niche_prop_mat) {
  
  niche_prop_mat <- acomp(niche_prop_mat) %>% as.matrix()
  
  # Generate ILR transformation
  baseILR <- ilrBase(x = niche_prop_mat,
                     method = "basic")
  cell_ilr <- as.matrix(ilr(niche_prop_mat, baseILR))
  colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))
  
  return(cell_ilr)
}

# Estimate explained variance with PCA regression

estimate_explvar <- function(ILR_mat) {
  
  pcs <- prcomp(x = ILR_mat) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_region_id") %>%
    left_join(patient_info %>% dplyr::select(patient_region_id, patient_group))
  
  pc_summary <- pc_summary %>%
    pivot_longer(-c("patient_region_id", "patient_group")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(aov_res = map(data, function(dat) {
      aov(value ~ patient_group,data = dat) %>%
        broom::tidy()
    }))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x))
  
  pc_summary <- pc_summary %>% 
    dplyr::select(aov_res) %>%
    unnest() %>%
    dplyr::filter(term != "Residuals") %>%
    left_join(prop_var) %>%
    ungroup() %>%
    dplyr::mutate(p.adj = p.adjust(p.value)) %>%
    arrange(p.adj)
  
  return(pc_summary)
  
}

estimate_explvar_time <- function(ILR_mat, early_only = T) {
  
  if(early_only) {
    
    useful_samples <- patient_info %>%
      dplyr::filter(! major_labl %in% c("FZ", "CTRL")) %>%
      pull(patient_region_id) %>%
      unique()
    
    #Remember that you filter patients before
    useful_samples <- useful_samples[useful_samples %in% rownames(ILR_mat)]
    
    ILR_mat <- ILR_mat[useful_samples, ]
    
  }
  
  pcs <- prcomp(x = ILR_mat) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_region_id") %>%
    left_join(patient_info %>% dplyr::select(patient_region_id, patient_group, major_labl, time))
  
  pc_summary <- pc_summary %>%
    pivot_longer(-c("patient_region_id", "patient_group", "major_labl", "time")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(aov_res = map(data, function(dat) {
      lm(value ~ time,data = dat) %>%
        broom::tidy()
    }))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x))
  
  pc_summary <- pc_summary %>% 
    dplyr::select(aov_res) %>%
    unnest() %>%
    dplyr::filter(term == "time") %>%
    left_join(prop_var) %>%
    ungroup() %>%
    dplyr::mutate(p.adj = p.adjust(p.value)) %>%
    arrange(p.adj)
  
  return(pc_summary)
  
}

estimate_explvar_area <- function(ILR_mat, myogenic_enriched = F) {
  
  if(myogenic_enriched) {
    
    useful_samples <- patient_info %>%
      dplyr::filter(! major_labl %in% c("FZ", "IZ")) %>%
      pull(patient_region_id) %>%
      unique()
    
    #Remember that you filter patients before
    useful_samples <- useful_samples[useful_samples %in% rownames(ILR_mat)]
    
    ILR_mat <- ILR_mat[useful_samples, ]
    
  }
  
  pcs <- prcomp(x = ILR_mat) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_region_id") %>%
    left_join(patient_info %>% dplyr::select(patient_region_id, major_labl))
  
  pc_summary <- pc_summary %>%
    pivot_longer(-c("patient_region_id", "major_labl")) %>%
    group_by(name) %>%
    nest() %>%
    mutate(aov_res = map(data, function(dat) {
      aov(value ~ major_labl,data = dat) %>%
        broom::tidy()
    }))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x))
  
  pc_summary <- pc_summary %>% 
    dplyr::select(aov_res) %>%
    unnest() %>%
    dplyr::filter(term != "Residuals") %>%
    left_join(prop_var) %>%
    ungroup() %>%
    dplyr::mutate(p.adj = p.adjust(p.value)) %>%
    arrange(p.adj)
  
  return(pc_summary)
  
}


# Generate clustering

plot_clust <- function(ILR_mat, explvar_val) {
  
  gex_hclust <- eclust(ILR_mat, "hclust", k = 3)
  
  # Make color palette
  
  color_palette <- tibble(patient_region_id = gex_hclust$labels[gex_hclust$order]) %>%
    left_join(patient_info[,c("patient_group", "patient_region_id")] %>% unique()) %>%
    left_join(tibble(patient_group = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"),
                     col = c("red", "darkgreen", "blue")))
  
  
  return(fviz_dend(gex_hclust, 
                 rect = TRUE, 
                 label_cols = color_palette$col,
                 k_colors = rep("black",3),
                 main = paste0("Expl. var = ", explvar_val)))
  
}

# Generate PCA instead of clustering - to do

plot_PCA <- function(ILR_mat, 
                     explvar_val, 
                     early_only = F,
                     res_dir) {
  
  
  
  
  if(early_only) {
    
    useful_samples <- patient_info %>%
      dplyr::filter(major_labl %in% c("IZ")) %>%
      pull(patient_region_id) %>%
      unique()
    
    #Remember that you filter patients before
    useful_samples <- useful_samples[useful_samples %in% rownames(ILR_mat)]
    
    ILR_mat <- ILR_mat[useful_samples, ]
    
  }
  
  pcs <- prcomp(x = ILR_mat) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_region_id") %>%
    left_join(patient_info %>% dplyr::select(patient_region_id, patient_group, major_labl, time))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x)) %>%
    dplyr::filter(name %in% c("PC1", "PC2"))
  
  prop_var <- sum(prop_var$expl_var)
  
  if(early_only) {
    
    plt <- ggplot(pc_summary, 
           aes(x = PC1, y = PC2, color = patient_group, label = time)) +
      geom_point() +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      theme(axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", 
                                        fill=NA, size=0.5)) +
      ggtitle(paste0("Prop. Expl. Var by time =", round(explvar_val, 2), "\n",
                     "Expl. Var by PCs = ", round(prop_var, 2)))
    
    pdf(paste0(res_dir, "/compPCA_early_only.pdf"), height = 6, width = 7)
    
    plot(plt)
    
    dev.off()
    
    write_csv(pc_summary, paste0(res_dir, "/compPCA_early_only.csv"))
    
    
  } else {
    
    plt <- ggplot(pc_summary, 
           aes(x = PC1, y = PC2, color = patient_group, label = patient_region_id)) +
      geom_point() +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      theme(axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", 
                                        fill=NA, size=0.5)) +
      ggtitle(paste0("Prop. Expl. Var by patient group =", round(explvar_val, 2), "\n",
                     "Expl. Var by PCs = ", round(prop_var, 2)))
    
    pdf(paste0(res_dir, "/compPCA.pdf"), height = 6, width = 7)
    
    plot(plt)
    
    dev.off()
    
    write_csv(pc_summary, paste0(res_dir, "/compPCA.csv"))

  }
 
  return(NULL)

}

plot_PCA_area <- function(ILR_mat, 
                     explvar_val, 
                     myogenic_enriched = F,
                     res_dir) {
  
  label <- "compPCA_area"
  
  if(myogenic_enriched) {
    
    useful_samples <- patient_info %>%
      dplyr::filter(! major_labl %in% c("FZ", "IZ")) %>%
      pull(patient_region_id) %>%
      unique()
    
    #Remember that you filter patients before
    useful_samples <- useful_samples[useful_samples %in% rownames(ILR_mat)]
    
    ILR_mat <- ILR_mat[useful_samples, ]
    
    label <- "compPCA_area_me"
    
  }
  
  pcs <- prcomp(x = ILR_mat) 
  
  pc_summary <- pcs$x %>%
    as.data.frame() %>%
    rownames_to_column("patient_region_id") %>%
    left_join(patient_info %>% dplyr::select(patient_region_id, patient_group, major_labl, time))
  
  prop_var <- tibble(expl_var = pcs$sdev/sum(pcs$sdev),
                     name = colnames(pcs$x)) %>%
    dplyr::filter(name %in% c("PC1", "PC2"))
  
  prop_var <- sum(prop_var$expl_var)

  plt <- ggplot(pc_summary, 
                  aes(x = PC1, y = PC2, 
                      color = major_labl, 
                      label = patient_region_id)) +
      geom_point() +
      ggrepel::geom_text_repel() +
      theme_minimal() +
      theme(axis.text = element_text(size = 10),
            panel.border = element_rect(colour = "black", 
                                        fill=NA, size=0.5)) +
      ggtitle(paste0("Prop. Expl. Var by area =", round(explvar_val, 2), "\n",
                     "Expl. Var by PCs = ", round(prop_var, 2)))
    
    pdf(paste0(res_dir, "/", label, ".pdf"), height = 6, width = 7)
    
    plot(plt)
    
    dev.off()
    
    write_csv(pc_summary, paste0(res_dir, "/", label, ".csv"))
  
  return(NULL)
  
}





# 4. Characterize with  cell-states

run_statecharacterization <- function(state_description, 
                                      niche_info,
                                      res_dir) {
  
  niche_description <- state_description %>%
    left_join(niche_info %>%
                dplyr::select(mol_niche, spot_id)) %>%
    na.omit() %>%
    select(-spot_id) %>%
    group_by(name) %>%
    nest() %>%
    mutate(wres = map(data, function(dat) {
      
      niches <- dat$mol_niche %>%
        unique() %>%
        set_names()
      
      map(niches, function(g) {
        
        test_data <- dat %>%
          mutate(test_group = ifelse(.data[["mol_niche"]] == g,
                                     "target", "rest")) %>%
          mutate(test_group = factor(test_group,
                                     levels = c("target", "rest")))
        
        wilcox.test(value ~ test_group, 
                    data = test_data,
                    alternative = "greater") %>%
          broom::tidy()
      }) %>% enframe("mol_niche") %>%
        unnest()
      
    }))
  
  wilcox_states <- niche_description %>%
    dplyr::select(wres) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
    dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
    dplyr::mutate(sign = ifelse(adj_pval < 0.05, "*", ""))
  
  write_csv(wilcox_states, file = paste0(res_dir, "/cell_type_enrich_wilcox.csv"))
  
  wilcox_states_plt <- wilcox_states %>%
    dplyr::mutate(log_adj_pval = ifelse(is.infinite(log_adj_pval), NaN, log_adj_pval)) %>%
    dplyr::mutate(log_adj_pval = ifelse(is.nan(log_adj_pval), max(log_adj_pval,na.rm = T), log_adj_pval))
  
  mapping_q <- quantile(wilcox_states_plt$log_adj_pval, 0.99)
  
  wilcox_states_plt <- wilcox_states_plt %>%
    mutate(log_adj_pval = ifelse(log_adj_pval >= mapping_q, mapping_q, log_adj_pval))
  
  niche_cs_plt <- ggplot(wilcox_states_plt, aes(x = name, y = ct_niche, fill = log_adj_pval)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_fill_gradient(na.value = "black",low = 'black',high = "yellow") +
    ylab("") +
    xlab("") +
    coord_equal()
  
  pdf(paste0(res_dir, "/cell_state_enrich_wilcox.pdf"), height = 6, width = 8)
  
  plot(niche_cs_plt)
  
  dev.off()
  
  return(NULL)
}

# 5. Characterize with cell-types

run_typecharacterization <- function(ct_description, 
                                     niche_info,
                                     res_dir) {
  
  cellprops_info <-  ct_description %>%
    left_join(niche_info %>%
                dplyr::select(mol_niche, spot_id, patient_region_id)) %>%
    na.omit()
  
  # Get median values per patient
  
  cell_props_summary_mol_pat <- cellprops_info %>%
    group_by(name, patient_region_id, mol_niche) %>%
    summarize(median_mol = median(value))
  
  # Check niches that are unique to some patients
  
  niche_selection <- cell_props_summary_mol_pat %>% 
    ungroup() %>% 
    group_by(mol_niche) %>%
    dplyr::select(patient_region_id, mol_niche) %>%
    unique() %>%
    summarise(n_samples = n()) %>%
    dplyr::filter(n_samples >= 5) %>%
    pull(mol_niche)
  
  cell_props_summary_mol <- cell_props_summary_mol_pat %>%
    dplyr::filter(mol_niche %in% niche_selection) %>%
    ungroup() %>%
    group_by(name, mol_niche) %>%
    summarize(median_mol_niche = median(median_mol))
  
  niche_summary_mat <- cell_props_summary_mol %>%
    pivot_wider(values_from = median_mol_niche, 
                names_from =  name, values_fill = 0) %>%
    column_to_rownames("mol_niche") %>%
    as.matrix()
  
  niche_order <- hclust(dist(niche_summary_mat))
  niche_order <- niche_order$labels[niche_order$order]
  
  ct_order <- hclust(dist(t(niche_summary_mat)))
  ct_order <- ct_order$labels[ct_order$order]
  
  # Find characteristic cell types of each niche
  # We have per patient the proportion of each cell-type in each niche
  
  run_wilcox_up <- function(prop_data) {
    
    prop_data_group <- prop_data[["mol_niche"]] %>%
      unique() %>%
      set_names()
    
    map(prop_data_group, function(g) {
      
      test_data <- prop_data %>%
        mutate(test_group = ifelse(mol_niche == g,
                                   "target", "rest")) %>%
        mutate(test_group = factor(test_group,
                                   levels = c("target", "rest")))
      
      wilcox.test(median_mol ~ test_group, 
                  data = test_data,
                  alternative = "greater") %>%
        broom::tidy()
    }) %>% enframe("mol_niche") %>%
      unnest()
    
  }
  
  wilcoxon_res <- cell_props_summary_mol_pat %>%
    dplyr::filter(mol_niche %in% niche_selection) %>%
    ungroup() %>%
    group_by(name) %>%
    nest() %>%
    mutate(wres = map(data, run_wilcox_up)) %>%
    dplyr::select(wres) %>%
    unnest() %>%
    ungroup() %>%
    mutate(p_corr = p.adjust(p.value))
  
  mean_ct_prop_plt <- cell_props_summary_mol %>%
    left_join(wilcoxon_res, by = c("mol_niche", "name")) %>%
    mutate(significant = ifelse(p_corr <= 0.15, "*", "")) %>%
    mutate(cell_type = factor(name, levels = ct_order),
           niche_mol = factor(mol_niche, levels = niche_order)) %>%
    ungroup() %>%
    group_by(cell_type) %>%
    mutate(scaled_pat_median = (median_mol_niche - mean(median_mol_niche))/sd(median_mol_niche)) %>%
    ungroup() %>%
    ggplot(aes(x = name, y = niche_mol, fill = scaled_pat_median)) +
    geom_tile() +
    geom_text(aes(label = significant)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
          legend.position = "bottom",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.y = element_text(size=12)) +
    scale_fill_gradient(high = "#ffd89b", low = "#19547b") +
    coord_equal() +
    ylab("") +
    xlab("")
  
  pdf(paste0(res_dir, "/cell_type_enrich_wilcox.pdf"), height = 5, width = 5)
  
  plot(mean_ct_prop_plt)
  
  dev.off()
  
  cell_props_summary_mol %>%
    left_join(wilcoxon_res, by = c("mol_niche", "name")) %>%
    mutate(significant = ifelse(p_corr <= 0.15, "*", "")) %>%
    mutate(cell_type = factor(name, levels = ct_order),
           niche_mol = factor(mol_niche, levels = niche_order)) %>%
    ungroup() %>%
    group_by(cell_type) %>%
    mutate(scaled_pat_median = (median_mol_niche - mean(median_mol_niche))/sd(median_mol_niche)) %>%
    ungroup() %>%
    write_csv(., file = paste0(res_dir, "/cell_type_enrich_wilcox.csv"))
  
  return(NULL)
  
}


run_characterization <- function(ct_description, 
                                 niche_info,
                                 res_dir,
                                 alias,
                                 height,
                                 width) {
  
  ct_description <-  ct_description %>%
    left_join(niche_info %>%
                dplyr::select(mol_niche, spot_id, patient_region_id)) %>%
    na.omit() %>%
    select(-spot_id) %>%
    group_by(name) %>%
    nest() %>%
    mutate(wres = map(data, function(dat) {
      
      niches <- dat$mol_niche %>%
        unique() %>%
        set_names()
      
      map(niches, function(g) {
        
        test_data <- dat %>%
          mutate(test_group = ifelse(.data[["mol_niche"]] == g,
                                     "target", "rest")) %>%
          mutate(test_group = factor(test_group,
                                     levels = c("target", "rest")))
        
        wilcox.test(value ~ test_group, 
                    data = test_data,
                    alternative = "greater") %>%
          broom::tidy()
      }) %>% enframe("mol_niche") %>%
        unnest()
      
    }))
  
  wilcox_types <- ct_description %>%
    dplyr::select(wres) %>%
    unnest() %>%
    ungroup() %>%
    dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
    dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
    dplyr::mutate(sign = ifelse(adj_pval < 0.005, "*", ""))
  
  ct_median_desc <- ct_description %>%
    dplyr::select(data) %>%
    unnest() %>%
    group_by(name, mol_niche) %>%
    summarise(median_prop = median(value)) %>%
    mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop))
  
  # Bind both dataframes
  
  niche_car_df <- left_join(wilcox_types, ct_median_desc) %>% na.omit()
  
  write_csv(niche_car_df, file = paste0(res_dir, "/", alias, ".csv"))
  
  # Order based on clustering of both rows and columns
  ct_median_desc_mat <- niche_car_df %>% dplyr::select(name, mol_niche, scaled_median_prop) %>%
    pivot_wider(names_from = mol_niche, values_from = scaled_median_prop) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  ct_sign_desc_mat <- niche_car_df %>% dplyr::select(name, mol_niche, sign) %>%
    pivot_wider(names_from = mol_niche, values_from = sign) %>%
    column_to_rownames("name") %>%
    as.matrix()
  
  
  niche_car_plt <- Heatmap(ct_median_desc_mat, name = "scaled comp", rect_gp = gpar(col = "black", lwd = 1),
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.text(sprintf(ct_sign_desc_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                           })
  
  pdf(paste0(res_dir, "/", alias, ".pdf"), height = height, width = width)
  
  draw(niche_car_plt)
  
  dev.off()
  
  return(NULL)
  
}















