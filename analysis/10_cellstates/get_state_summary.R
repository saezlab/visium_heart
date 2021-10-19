# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

# Here we generalize the analysis of states from a single cell type object

library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(HDF5Array)
library(uwot)
library(lme4)
library(lmerTest)
library(progeny)
library(dorothea)
library(viper)
library(fgsea)
library(cowplot)
source("./analysis/utils/pseudobulk_utils.R")

# This analysis must be done at patient_region level to take into account repetitions

summarize_state <- function(ct_folder, ct_alias, exception, 
                            perc_thrsh, n_samples_filt,
                            background_exclude = NULL){
  
  # Define outs -------------------------------------------------------
  gene_filter_out <- paste0(ct_folder, ct_alias, "/gene_filter_df.rds")
  pb_path <- paste0(ct_folder, ct_alias, "/ps_", ct_alias, "_states.rds")
  state_prop_results <- paste0(ct_folder, ct_alias, "/state_proportion_kwtest.txt")
  mlm_res <-  paste0(ct_folder, ct_alias, "/state_mrks.txt")
  state_gene_list <- paste0(ct_folder, ct_alias, "/state_genelist.rds")
  gsea_out <- paste0(ct_folder, ct_alias, "/state_gsea.txt")
  progeny_out <- paste0(ct_folder, ct_alias, "/state_progeny.txt")
  de_out <- paste0(ct_folder, ct_alias, "/state_de.txt")
  tf_out <- paste0(ct_folder, ct_alias, "/state_tf.txt")
  all_out <- paste0(ct_folder, ct_alias, "/states_funcomics.pdf")
  cor_plt <- paste0(ct_folder, ct_alias, "/states_cor.pdf")
  reg_net <- paste0("./reg_nets/processed/", ct_alias, ".txt")
  
  # Read object and filter gene expression ---------------------------------------
  print("Reading scell object and filtering lowly expressed genes")
  sc_data_file <- paste0(ct_folder, ct_alias, "/", ct_alias, "_states_sce/")
  sc_data <- loadNfilter_sce(sc_data_file, gene_filter_out, exception)
  sc_meta <- colData(sc_data) %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    left_join(patient_annotations, 
              by = c("orig.ident" = "sample_id"))
  
  # Get UMAP for visualization
  sc_data_umap <- reducedDim(sc_data, "UMAP") %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    left_join(sc_meta, by = "cell_id") 
  
  # Get state specific genes ----------------------------------------
  niche_info <- readRDS(gene_filter_out)
  
  # Estimate state proportions while filtering for pat specific states
  
  cell_state_prop <- get_pat_state_props(sc_meta, 
                      n_samples_filt = n_samples_filt, 
                      min_pat_perc = 0.01)
  
  filtered_states <- cell_state_prop$opt_state %>% unique()
  
  # N cells after filtering
  n_cells_plot <- cell_state_prop %>%
    dplyr::select(patient_id, n_cells_sample) %>%
    unique() %>%
    ggplot(aes(x = n_cells_sample, y = patient_id)) +
    geom_bar(stat = "identity") +
    ylab("sample") +
    xlab("number of called cells")
  
  # Prop of states per sample
  sample_states_dots <- cell_state_prop %>%
    ggplot(aes(x = opt_state, y = patient_id, size = props_state)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
  
  # UMAP
  sc_data_umap_plt <- sc_data_umap %>%
    dplyr::filter(opt_state %in% filtered_states) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2, color = opt_state)) +
    geom_point(size = 0.2) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  # Compare compositions between patient groups ---------------------
  
  kwallis_res <- compare_state_comps(cell_state_prop,state_prop_results)
  
  kwallis_winners <- dplyr::filter(kwallis_res, p_corr < 0.15) %>%
    pull(opt_state)
  
  cell_state_prop_plt <- cell_state_prop %>%
    left_join(patient_annotations, by = "patient_id") %>%
    mutate(opt_state = ifelse(opt_state %in% kwallis_winners, 
                              paste0(opt_state, "*"), opt_state)) %>%
    ggplot(aes(x = patient_group, y = props_state, color = patient_group)) +
    geom_boxplot() +
    facet_wrap(.~opt_state, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # Pseudobulk analysis -----------------------------------------------
  
  # Get the data, and filter for states that we find useful
  
  pb_data <- get_pseudobulk_mat(pb_path = pb_path,
                     patient_annotations = patient_annotations,
                     filtered_states = filtered_states,
                     niche_info = niche_info)
  
  # Second, filter out genes that are background
  
  if(!is.null(background_exclude)) {
    
    marker_list_cts <- readRDS("./markers/pb_ct_marker_list.rds")
    
    exclude_cts <- names(marker_list_cts)
    
    exclude_cts <- exclude_cts[!(exclude_cts %in% background_exclude)]
    
    marker_list_cts <- marker_list_cts[exclude_cts] %>%
      unlist() %>% 
      unique()
    
    current_genes <- pb_data$gex %>%
      rownames()
    
    current_genes <- current_genes[!(current_genes %in% marker_list_cts)]
    
    pb_data$gex <- pb_data$gex[current_genes, ]
    
  }
  
  # Here we also need to exclude background markers
  
  pb_data <- filt_notinf_pbprofiles(pb_data = pb_data, pb_perc = 0.01)
  
  pb_meta <- pb_data$meta
  pb_gex <- pb_data$gex
  
  # Normalizing
  pb_gex_norm <- pb_gex %>%
    cpm_norm(expression_matrix = .)
  
  # Long object
  pb_cpm_long <- pb_gex_norm %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "col_id", values_to = "expr") %>%
    left_join(pb_meta, by = "col_id")
  
  # Run differential expression analysis -------------------------------------
  print("Running pseudobulk differential analysis")
  
  degs_ext <- run_edgeR(pb_gex = pb_gex,
                        pb_meta = pb_meta) %>%
    dplyr::rename("state" = name)
  
  degs_ext %>%
    dplyr::filter(logFC > 0) %>% 
  write.table(., col.names = T, 
              row.names = F, quote = F, sep = "\t",
              file = mlm_res)
  
  # Correlation of states
  cor_states <- degs_ext %>%
    dplyr::select(state, gene, logFC) %>%
    pivot_wider(names_from = state, 
                values_from = logFC,
                values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    cor(.,method = "spearman") 
  
  order_states <- hclust(as.dist(1-cor_states))
  order_states <- order_states$labels[order_states$order]
  
  cor_states <- cor_states %>%
    as.data.frame() %>%
    rownames_to_column("state_a") %>%
    pivot_longer(-state_a, names_to = "state_b", values_to = "spearman_cor") %>%
    mutate(state_a = factor(state_a, levels = order_states),
           state_b = factor(state_b, levels = order_states))
  
  cor_plt <- ggplot(cor_states, aes(x = state_b, y = state_a, fill = spearman_cor)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=90, 
                                     hjust = 1, 
                                     vjust = 0.5)) +
    scale_fill_gradient2(limits = c(-1,1))
  
  # Plotting marker genes
  
  plot_genes <- degs_ext %>% 
    mutate(corr_pval = FDR) %>%
    dplyr::filter(corr_pval < 0.15) %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(state, -logFC) %>%
    group_by(state) %>%
    dplyr::slice(1:5) %>%
    pull(gene) %>%
    unique()
  
  #Showing the expression of marker genes 
  pb_mrkr_plt <- pb_cpm_long %>% 
    dplyr::select(gene, col_id, expr, opt_state) %>%
    dplyr::filter(gene %in% plot_genes) %>%
    arrange(gene, opt_state) %>%
    group_by(gene) %>%
    mutate(expr = (expr - mean(expr)) / sd(expr),
           gene = factor(gene,
                         levels = plot_genes),
           col_id = factor(col_id,
                           levels = col_id)) %>%
    ggplot(aes(x = gene, y = col_id, fill = expr)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient2()
  
  pb_mrkr_plt_ts <- degs_ext %>% 
    dplyr::filter(gene %in% plot_genes) %>%
    mutate(gene = factor(gene,
                         levels = plot_genes)) %>%
    ggplot(aes(y = gene, x = state, fill = logFC)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_gradient2()
  
  # Get representative genes for later:
  # This is to enrich in visium slides
  
  state_genes <- degs_ext %>% 
    arrange(state, -logFC) %>%
    dplyr::filter(FDR <= 0.15, logFC > 0) %>%
    group_by(state) %>%
    dplyr::slice(1:200) %>%
    dplyr::select(gene, state) %>%
    nest() %>%
    deframe() %>%
    map(., ~.x[[1]])
  
  saveRDS(state_genes, state_gene_list)
  
  # Funcomics -------------------------------
  
  # GSEA -------
  gsea_res <- run_gsea(degs_ext = degs_ext,
                       gsea_out = gsea_out)
  
  # PROGENY -------
  progeny_res <- run_progeny_t(degs_ext = degs_ext, progeny_out = progeny_out)
  # order of pathways in plot
  progeny_order <- hclust(dist(t(progeny_res)))
  progeny_order <- progeny_order$labels[progeny_order$order]
  
  progeny_plt <- progeny_res %>%
    as.data.frame() %>%
    rownames_to_column("state") %>%
    pivot_longer(-state, names_to = "pathway", values_to = "progeny_scores") %>%
    mutate(pathway = factor(pathway, levels = progeny_order)) %>%
    ggplot(aes(x = pathway, y = state, fill = progeny_scores)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # Dorothea + wmean -------
  
  dorothea_res <- run_dorothea_t(degs_ext = degs_ext,
                                 tf_out = tf_out,
                                 reg_net = reg_net)
  
  dorothea_res_sign <- dorothea_res %>%
    arrange(condition, -score) %>%
    group_by(condition) %>%
    dplyr::slice(1:10) 
  
  dorothea_plt <- dorothea_res %>%
    dplyr::select(condition, source, score) %>%
    dplyr::filter(source %in% dorothea_res_sign$source) %>%
    dplyr::mutate(source = factor(source,
                                  levels = unique(dorothea_res_sign$source))) %>%
    ggplot(aes(y = source, x = condition, fill = score)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # Final panel
  
  upper_panel <- cowplot::plot_grid(n_cells_plot, sample_states_dots, 
                                    sc_data_umap_plt, ncol = 3,
                                    rel_widths = c(0.2, 0.4, 0.5), align = "hv")
  
  lower_panel_left <- cowplot::plot_grid(cell_state_prop_plt, progeny_plt, ncol = 1, rel_heights = c(0.5,0.4))
  
  lower_panel_right <- cowplot::plot_grid(pb_mrkr_plt_ts, dorothea_plt, nrow = 1)
  
  lower_panel <- cowplot::plot_grid(lower_panel_left, lower_panel_right, nrow = 1, rel_widths = c(0.5, 0.4))
  
  pdf(all_out, height = 18, width = 18)
  all_panels <- (cowplot::plot_grid(upper_panel,lower_panel, ncol = 1, rel_heights = c(0.35,0.5)))
  plot(all_panels)
  plot(cor_plt)
  dev.off()
  
  return(all_panels)
  
}

# This function filters lowly expressed genes based on min number of cells in state of interest
# This also saves a list of state specific genes in gene_filter_out
# returns the sce object
loadNfilter_sce <- function(sc_data_file, gene_filter_out, exception) {
  # Load
  sc_data <- loadHDF5SummarizedExperiment(sc_data_file)
  # Add label to state
  sc_data[["opt_state"]]  <- paste0("state", sc_data[["opt_state"]])
  
  print("Excluding exceptions")
  
  # Filter niches that are trash -------------------------------
  if(is.na(exception)) {
    
    sc_meta <- colData(sc_data) %>%
      as.data.frame() %>%
      rownames_to_column("cell_id")
    
  } else {
    
    sc_meta <- colData(sc_data) %>%
      as.data.frame() %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(! opt_state %in% paste0("state", exception))
    
    sc_data <- sc_data[, sc_meta[["cell_id"]]]
  }
  
  # Calculate the number of cells per state
  # The minimum will represent the filtering of lowly expressed genes
  
  print("Excluding lowly expressed genes")
  
  n_cell_df <- colData(sc_data) %>% 
    as.data.frame() %>%
    group_by_at("opt_state") %>%
    summarise(n_spots = n())
  
  min_spots <- n_cell_df %>%
    pull(n_spots) %>%
    min()
  
  # First let's filter all genes that are extremely lowly expressed
  expressed_genes <- counts(sc_data) > 0
  gene_ix <- rowSums(expressed_genes) > min_spots
  sc_data <- sc_data[gene_ix, ]
  expressed_genes <- expressed_genes[gene_ix, ]
    
  print("Keeping genes expressed in at least perc_thrsh in a state")
  
  # Then for each niche we will test specifically if the gene can be considered expressed
  # Get cell_ids that belong to a class
  niche_info <- colData(sc_data) %>% 
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    dplyr::select_at(c("cell_id", "opt_state")) %>%
    group_by_at("opt_state") %>%
    nest() %>%
    dplyr::rename("cell_ids" = data) %>%
    dplyr::mutate(cell_ids = map(cell_ids, ~ .x[[1]])) %>%
    left_join(n_cell_df) %>%
    dplyr::mutate(min_spots = (n_spots * perc_thrsh) %>% floor()) %>%
    dplyr::mutate(selected_genes = map2(cell_ids, min_spots, function(cids, mspots) {
      gene_ix <- rowSums(expressed_genes[, cids]) > mspots
      return(names(gene_ix[gene_ix]))
    }))
  
  # This data frame contains the genes that can be considered as expressed in a niche
  # under our thresholds
  saveRDS(niche_info, gene_filter_out)
  
  return(sc_data)
}

# This function filters states that aren't representative of all samples

get_pat_state_props <- function(sc_meta, n_samples_filt, min_pat_perc = 0.01) {
  
  print("Calculating proportions of states per patient")
  
  cell_state_prop <- sc_meta %>%
    dplyr::group_by_at(c("patient_id", "opt_state")) %>%
    summarise(n_cells_state = n()) %>%
    ungroup() %>%
    group_by(patient_id) %>%
    mutate(n_cells_sample = sum(n_cells_state)) %>%
    ungroup() %>%
    mutate(props_state = n_cells_state/n_cells_sample)
  
  print("Calculating states that represent at least min_pat_perc ")
  
  state_counts <- cell_state_prop %>%
    dplyr::filter(props_state > min_pat_perc) %>% # Statee represents at least 1% of your pop?
    group_by_at("opt_state") %>%
    summarize(n_samples = n())
  
  print("Keeping states represented in at least n_samples_filt, controlled by min_pat_perc")
  
  filtered_states <- state_counts %>%
    dplyr::filter(n_samples >= n_samples_filt) %>%
    pull("opt_state")
  
  cell_state_prop <- cell_state_prop %>%
    filter(get("opt_state") %in% filtered_states)
  
  return(cell_state_prop)
}

# Perform Kruskall-Wallis to identify changes in composition
# Distributions

compare_state_comps <- function(cell_state_prop, state_prop_results) {
  
  print("compare proportions between patient groups")
  
  kwallis_res <- cell_state_prop %>%
    left_join(patient_annotations, by = "patient_id") %>%
    group_by_at("opt_state") %>%
    nest() %>%
    mutate(kwallis_test = map(data, function(dat){
      broom::tidy(kruskal.test(props_state ~ patient_group, data = dat))
    })) %>%
    dplyr::select(kwallis_test) %>%
    unnest() %>%
    ungroup() %>%
    mutate(p_corr = p.adjust(p.value))
  
  kwallis_res %>% 
    write.table(file = state_prop_results, quote = F, 
                col.names = T, row.names = F, sep = "\t")
  
  return(kwallis_res)
  
}

# Here we generate the pseudobulk profile per patient rather than sample
# We used filtered states and genes that we know we should evaluate
get_pseudobulk_mat <- function(pb_path, patient_annotations, 
                               filtered_states, niche_info) {
  
  print("getting pseudobulk matrix and filtering with genes and niches of interest")
  
  pb_data <- readRDS(pb_path)[[1]][["gex"]]
  
  pb_meta <- colData(pb_data) %>% 
    as.data.frame() %>% 
    left_join(patient_annotations, 
              by = c("orig.ident" = "sample_id"))
  
  patient_cells <- pb_meta %>%
    group_by(patient_id, opt_state) %>%
    summarize(ncells = sum(ncells))
  
  # This summarizes the info into patients
  pb_data <- sumCountsAcrossCells(assay(pb_data), 
                                  DataFrame(pb_meta[, c("patient_id", "opt_state")]))
  
  pb_meta <- colData(pb_data) %>% 
    as.data.frame() %>% 
    dplyr::select(-ncells) %>%
    left_join(patient_cells, by = c("patient_id","opt_state")) %>%
    left_join(patient_annotations %>%
                dplyr::select(patient_group, region_novel,patient_id) %>%
                unique(), 
              by = c("patient_id"))
  
  pb_meta[, "opt_state"] <- paste0("state", pb_meta[, "opt_state"])
  pb_meta[, "col_id"] <- paste0(pb_meta$patient_id, "_", pb_meta$opt_state)
  pb_gex <- assay(pb_data)
  colnames(pb_gex) <- pb_meta[, "col_id"]
  
  # Get states that we know are interesting and genes that we should evaluate
  niche_filter_ix <- pb_meta[, "opt_state"] %in% filtered_states
  pb_meta <- pb_meta[niche_filter_ix, ]
  pb_gex <- pb_gex[niche_info$selected_genes %>% unlist() %>% unique(),
                   niche_filter_ix]
  
  return(list("gex" = pb_gex, "meta" = pb_meta))
  
}

# This filters profiles coming from a meaningless population
# of cells of a single sample
# This could also be done at the number of cells...

filt_notinf_pbprofiles <- function(pb_data, pb_perc = 0.01, by = "ncells"){
  
  print("Excluding profiles that represent less than pb_perc in a patient")
  
  pb_gex <- pb_data$gex
  
  pb_meta <- pb_data$meta
  
  if(by == "ncells"){
    
    pb_meta$max_counts <- colMaxs(pb_gex)
    pb_meta$keep <- ifelse(pb_meta$ncells < 50 &
                             pb_meta$max_counts < 1000, 
                             FALSE, TRUE)
    t_ix <- which(pb_meta$keep == TRUE)
    pb_meta <- pb_meta[t_ix, ]
    pb_gex <- pb_gex[, t_ix]
    
    return(list("gex" = pb_gex, "meta" = pb_meta))
    
  } else{
    
    # Then rejecting all pseudobulk profiles that come from less than certain percent of population
    pb_meta_props <- pb_meta %>%
      group_by(patient_id) %>%
      mutate(ncells_pat = sum(ncells)) %>%
      mutate(propcells_pat = ncells/ncells_pat)
    
    perc_filt <- which(pb_meta_props$propcells_pat >= pb_perc)
    pb_gex <- pb_gex[, perc_filt]
    pb_meta <- pb_meta[perc_filt, ]
    
    return(list("gex" = pb_gex, "meta" = pb_meta))
  }
  
}

run_edgeR <- function(pb_gex, pb_meta) {
  
  cts <- set_names(pb_meta$opt_state %>% unique)
  
  de_res <- map(cts, function(ct) {
    print(ct)
    
    ct_meta_data <- pb_meta %>%
      mutate(test_column = ifelse(opt_state == ct, ct, "rest"))
    
    dat <- DGEList(pb_gex, samples = DataFrame(ct_meta_data))
    
    keep <- filterByExpr(dat, group = ct_meta_data$test_column)
    
    dat <- dat[keep,]
    
    dat <- calcNormFactors(dat)
    
    design <- model.matrix(~factor(test_column,
                                   levels = c("rest",ct)), dat$samples)
    
    colnames(design) <- c("int", ct)
    
    dat <- estimateDisp(dat, design)
    
    fit <- glmQLFit(dat, design, robust=TRUE)
    
    res <- glmQLFTest(fit, coef=ncol(design))
    
    de_res <- topTags(res, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
    
    return(de_res)
    
  })
  
  de_res <- de_res %>% 
    enframe() %>%
    unnest()
  
  return(de_res)
  
}

# Here we run linear mixed models to identify marker genes

run_mlm <- function(pb_cpm_long, mlm_res) {
  print("Running differential expression analysis with mlm")
  # You run this once
  pb_cpm_long <- pb_cpm_long %>%
    group_by(gene) %>%
    nest() %>%
    mutate(mxeff_models = map(data, function(x) {
      states <- x %>% pull(opt_state) %>% unique() %>% set_names()
      state_column <- x %>% pull(opt_state)

      map(states, function(niche_id) {

        niche_data <- x %>%
          dplyr::mutate(niche = ifelse(opt_state == niche_id,
                                       "test_niche",
                                       "reference_niche")) %>%
          dplyr::mutate(niche = factor(niche,
                                       levels = c("reference_niche",
                                                  "test_niche")))

        fit <- lmerTest::lmer(expr ~ niche + (1 | patient_id), data = niche_data)

      }) %>% enframe()

    }))

  test_resuls <- pb_cpm_long %>%
    dplyr::select(gene, mxeff_models) %>%
    unnest() %>%
    mutate(model_coefficients = map(value, function(x) {
      mem_res <- summary(x)
      mem_res$coefficients  %>%
        as.data.frame() %>%
        rownames_to_column("term") %>%
        pivot_longer(-term)
    })) %>%
    mutate(anova_pval = map(value, function(x) {
      aov_res <- anova(x)
      aov_res$`Pr(>F)`
    }))

  degs <- test_resuls %>%
    dplyr::select(-value) %>%
    unnest(anova_pval) %>%
    dplyr::select(-model_coefficients)

  degs_ext <- test_resuls %>%
    dplyr::select(-c("value", "anova_pval")) %>%
    dplyr::rename("state" = name) %>%
    unnest(model_coefficients) %>%
    dplyr::filter(term == "nichetest_niche",
                  name %in% c("t value", "Pr(>|t|)")) %>%
    dplyr::select(-term) %>%
    mutate(name = ifelse(name == "t value",
                         "t_value",
                         "pval")) %>%
    pivot_wider(names_from = name,
                values_from = value) %>%
    ungroup() %>%
    arrange(state,-t_value)

  write.table(degs_ext, col.names = T, 
              row.names = F, quote = F, sep = "\t",
              file = mlm_res)
  
  return(degs_ext)
  
}

# Funcomics functions

run_gsea <- function(degs_ext, gsea_out) {
  print("running gsea")
  
  gene_sets <- readRDS(file = "./markers/Genesets_Dec19.rds")[["MSIGDB_CANONICAL"]]
  
  gls <- dplyr::select(degs_ext, 
                       state, gene, 
                       FDR, logFC) %>%
    dplyr::mutate(pval = -log10(FDR)) %>%
    mutate(pval = ifelse(is.infinite(pval), NA, pval))
  
  max_pval <- max(gls$pval,na.rm = T)
  
  gls <- gls %>%
    dplyr::mutate(pval = ifelse(is.na(pval), max_pval, pval)) %>%
    dplyr::mutate(rank_ord = pval * logFC)
  
  # First GSEA
  gsea_res <- gls %>%
    dplyr::select(-c("pval", "FDR", "logFC")) %>%
    group_by(state) %>%
    nest() %>%
    mutate(data = map(data, function(x) {
      set_names(x$rank_ord, x$gene) %>%
        na.omit()
    })) %>%
    rename("gls_rank" = data) %>%
    mutate(gsea_stats = map(gls_rank, function(t_vals) {
      fgsea(pathways = gene_sets,stats = t_vals)
    }))
  
  gsea_res <- gsea_res %>%
    select(gsea_stats) %>%
    unnest() %>%
    dplyr::filter(padj < 0.15) %>%
    arrange(state,-(NES))
  
  gsea_res %>% 
    dplyr::select(-leadingEdge) %>%
    write.table(., sep = "\t",
                row.names = F, col.names = T,
                quote = F, file = gsea_out)
  
  return(gsea_res)
}

run_progeny_t <- function(degs_ext, progeny_out) {
  print("running progeny")
  
  gls <- dplyr::select(degs_ext, 
                       state, gene, 
                       FDR, logFC) %>%
    dplyr::mutate(pval = -log10(FDR)) %>%
    mutate(pval = ifelse(is.infinite(pval), NA, pval))
  
  max_pval <- max(gls$pval,na.rm = T)
  
  gls <- gls %>%
    dplyr::mutate(pval = ifelse(is.na(pval), max_pval, pval)) %>%
    dplyr::mutate(rank_ord = pval * logFC)
  
  progeny_res <- gls %>% 
    dplyr::select(-c("pval", "FDR", "logFC")) %>%
    pivot_wider(names_from = state,
                values_from = rank_ord,
                values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    progeny(.,scale = T,top = 500)
  
  progeny_res %>%
    as.data.frame() %>%
    rownames_to_column("state") %>%
    pivot_longer(-state, names_to = "pathway", values_to = "progeny_scores") %>%
    write.table(., sep = "\t",
                row.names = F, col.names = T,
                quote = F, file = progeny_out)
  
  return(progeny_res)
}

run_dorothea_t <- function(degs_ext, tf_out, reg_net) {
  print("running tf act inference")
  
  min_targets <- 5
  
  regulons <- read_table2(reg_net)
  
  #data(dorothea_hs, package = "dorothea")
  
  #regulons <- dorothea_hs %>% 
  #  dplyr::filter(confidence != "E")
  
  regulons <- regulons %>% 
    dplyr::filter(target %in% (degs_ext$gene %>% unique()))
  
  filtered_regulons <- regulons %>% 
    dplyr::select(source, target) %>%
    group_by(source) %>%
    nest() %>%
    mutate(data = map(data, ~ .x[[1]])) %>%
    mutate(n_targets = map_dbl(data, length)) %>%
    dplyr::filter(n_targets >= 5) %>%
    dplyr::select(-n_targets) %>%
    deframe() %>%
    names()
  
  regulons <- regulons %>%
    dplyr::filter(source %in% filtered_regulons)
  
  # Generate matrix of t-values
  
  gls <- dplyr::select(degs_ext, 
                       state, gene, 
                       FDR, logFC) %>%
    dplyr::mutate(pval = -log10(FDR)) %>%
    mutate(pval = ifelse(is.infinite(pval), NA, pval))
  
  max_pval <- max(gls$pval,na.rm = T)
  
  gls <- gls %>%
    dplyr::mutate(pval = ifelse(is.na(pval), max_pval, pval)) %>%
    dplyr::mutate(rank_ord = pval * logFC)
  
  t_mat <-  gls %>% 
    dplyr::select(-c("pval", "FDR", "logFC")) %>%
    pivot_wider(names_from = state,
                values_from = rank_ord,
                values_fill = 0) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  dorothea_res <- decoupleR::run_wmean(mat = t_mat, 
                       network = regulons,
                       .source = "source",
                       .target = "target",
                       .mor = "mor",
                       .likelihood = "likelihood",
                       times = 10000)
  
  dorothea_res <- dorothea_res %>%
    dplyr::filter(statistic == "norm_wmean") %>%
    group_by(condition) %>%
    mutate(padj = p.adjust(p_value))
  
  dorothea_res %>% arrange(condition, -abs(score)) %>%
    write.table(., sep = "\t",
                row.names = F, col.names = T,
                quote = F, file = tf_out)
  
  return(dorothea_res)
  
}
