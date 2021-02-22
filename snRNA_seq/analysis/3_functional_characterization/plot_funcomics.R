# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we plot all the funcomics results
library(optparse)
library(tidyverse)
library(Seurat)
source("./utils/dea.R")

# Argument definition ---------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "scell seurat file"),
  make_option(c("--dea_data_path"), 
              action ="store", 
              default = NULL, 
              type = 'character',
              help = "dea object obtained from estimate_dea"),
  make_option(c("--group_class"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "Identity to characterize"),
  make_option(c("--nfeats"), 
              action ="store", 
              default = 10, 
              type = 'double',
              help = "number of top features to show per cluster"),
  make_option(c("--ngenes_ORA"), 
              action ="store", 
              default = 10, 
              type = 'double',
              help = "minimum number of genes for ORA"),
  make_option(c("--pvalue_ORA"), 
              action ="store", 
              default = 0.15, 
              type = 'double',
              help = "minimum number of genes for ORA"),
  make_option(c("--test_assays"), 
              action= "store", 
              default = NULL, 
              type = 'character',
              help = "assays to test separated by commas"),
  make_option(c("--pvalue"), 
              action= "store", 
              default = 0.001, 
              type = 'double',
              help = "log fold change threshold"),
  make_option(c("--lfc"), 
              action= "store", 
              default = 0.5, 
              type = 'double',
              help = "log fold change threshold"),
  make_option(c("--out_path"), 
              action= "store", 
              default = "default", 
              type = 'character',
              help = "path prefix where data will be generated path/id")
)

# Parse the parameters ---------------------------------------------------------------------------------
opt <- parse_args(OptionParser(option_list = option_list))

cat("[INFO] Input parameters\n", file = stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Read data -------------------------------------------------------------------------------------------
integrated_data <- readRDS(data_path)
Idents(integrated_data) <- group_class

dea_res <- readRDS(dea_data_path)
# This is hardcoded somewhere
gsets <- readRDS("./markers/Genesets_Dec19.rds")
gsets <- gsets$MSIGDB_CANONICAL
# Get assays to test -------------------------------------------------------------------------------------------
test_assays <- unlist(strsplit(test_assays, ","))

# Create df to analyze
dea_res_df <- dea_res %>%
  enframe(name = "assay",
          value = "dea") %>%
  dplyr::filter(assay %in% test_assays) %>%
  dplyr::mutate(out_dat = paste0(out_path, "_", assay, ".txt"),
                out_fig = paste0(out_path, "_", assay, ".pdf"))

# Function to filter dea_res
filter_def <- function(dea_df, out_dat_path) { 
  
  filtered_df <- dplyr::filter(dea_df,
                            p_val_adj <= pvalue,
                            avg_logFC > lfc) %>%
    arrange(cluster,-avg_logFC)
  
  write.table(x = filtered_df,
              file = out_dat_path,
              quote = F,
              row.names = F,
              col.names = T, 
              sep = "\t")
  
  return(filtered_df)
  
}

# Function to print dotplots
plot_dots <- function(assay, dea_df) { 
  
  if(assay == "RNA") {
    plt_df <- dea_df %>%
      dplyr::filter(!grepl("MT-",gene)) %>%
      group_by(cluster) %>%
      dplyr::slice(1:nfeats)
  } else {
    plt_df <- dea_df %>%
      group_by(cluster) %>%
      dplyr::slice(1:nfeats)
  }
  
  feat_list <- unique(plt_df$gene)
  
  gene_dots <- DotPlot(integrated_data,
                       features = feat_list,
                       assay = assay) +
    coord_flip() + theme(axis.text.y = element_text(size=7),
                         axis.text.x = element_text(size=8,angle = 90,
                                                    vjust = 0.5),
                         axis.title = element_blank(),
                         legend.text = element_text(size=8),
                         legend.title = element_text(size=8))
}

# Main ---------------------------------------------
dea_res_df <- dea_res_df %>%
  mutate(dea = map2(dea, out_dat, filter_def)) %>%
  mutate(dplots = map2(assay, dea, plot_dots))

# Plot dots ---------------------------------
walk2(dea_res_df$dplots, dea_res_df$out_fig, function(x, y) { 
  
  pdf(file = y, height = 9, width = 6)
  
  print(x)
  
  dev.off()

})

# Generate ora -----------------------------

if("RNA" %in% dea_res_df$assay) {
  
  genes_df <- dea_res_df %>%
    dplyr::filter(assay == "RNA") %>%
    dplyr::select(dea) %>% unnest()
  
  genes <- genes_df %>%
    dplyr::select(cluster, gene) %>%
    group_by(cluster) %>%
    mutate(cluster = paste0("state_",as.character(cluster))) %>%
    nest() %>%
    deframe()
  
  genes <- map(genes, ~ .x[[1]])
  
  # filter mitochondrial genes
  
  genes <- map(genes, function(x) { 
    
    x[!grepl("MT-", x)]
    
    })
  
  # states to test for ORA
  states <- names(genes)
  states <- states[map_int(genes, length) > ngenes_ORA]
  genes <- genes[states]
  
  # Run ORA
  ora_res <- map(genes, GSE_analysis, Annotation_DB = gsets)

  ora_res <- map(ora_res, function(x) { 
    x %>% 
      dplyr::filter(corr_p_value <= pvalue_ORA)
    }) %>% 
    enframe(name = "state") %>%
    unnest()
  
  # Write ORA
  write.table(x = ora_res,
              file = paste0(out_path, "_ORA", ".txt"),
              quote = F,
              row.names = F,
              col.names = T, 
              sep = "\t")
  
  # Heatmap of ORA
  
  short_ora_res <- ora_res %>%
    group_by(state) %>%
    dplyr::slice(1:15) %>%
    ungroup()
  
  hmapora <- short_ora_res %>%
    ggplot(aes(x = state,
               y = factor(gset, 
                          levels = unique(short_ora_res$gset)), 
               fill = -log10(corr_p_value))) + 
    geom_tile(na.rm = T) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90,vjust=1),
          axis.text.y = element_text(hjust = 1)) +
    scale_fill_gradient(
      low = "black",
      high = "yellow",
      limits = c(0,15)) +
    ylab("")
  
  pdf(file = paste0(out_path, "_ORA", ".pdf"), 
      height = 9, width = 12)
  
  print(hmapora)
  
  dev.off()
    
}
