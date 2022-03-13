# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate funcomics
#' of niches from pseudobulk profiles

library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)
library(progeny)

source("./analysis/utils/pseudobulk_utils.R")

pseudobulk_data <- readRDS("./processed_visium/integration/pb_nichepat_dat.rds")
patient_meta <- readRDS("./markers/visium_patient_anns_revisions.rds")
sc_data_file <- readRDS()

counts <- assay(pseudobulk_data)
meta_data <- colData(pseudobulk_data) %>%
  as.data.frame()

meta_data$max_counts <- colMaxs(counts)

# We will exclude all profiles that are coming from less than 10 cells and a min of 1000 max counts
meta_data$keep <- ifelse((meta_data$ncells > 10) &
                           (meta_data$max_counts > 1000), 
                         TRUE, FALSE)

t_ix <- which(meta_data$keep == TRUE)

meta_data <- meta_data[t_ix, ]
counts <- counts[, t_ix]

meta_data$col_id <- paste0(meta_data$patient_id,
                           "_",
                           meta_data$niche)

colnames(counts) <- meta_data$col_id

# Filter for highly expressed genes
niche_info <- readRDS("./results/niche_mapping/expressed_genes_niche.rds")
useful_genes <- unlist(niche_info$selected_genes ) %>% unique()
counts <- counts[useful_genes, ]

cts <- set_names(meta_data$niche %>% unique)

de_res <- map(cts, function(ct) {
  print(ct)
  ct_meta_data <- meta_data %>%
    mutate(test_column = ifelse(niche == ct, ct, "rest"))
  
  dat <- DGEList(counts, samples = DataFrame(ct_meta_data))
  
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

de_res %>%
  dplyr::filter(logFC > 0) %>%
  arrange(name, FDR, - logFC) %>%
  write_csv(file = "./results/niche_mapping/edgeR_nichemrkrs.csv")

# Get top 10 genes per cell
pos_genes <- de_res %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  dplyr::filter(PValue <= 0.25) %>%
  arrange(name, - (logFC)) %>%
  group_by(name) %>%
  dplyr::slice(1:10) %>%
  pull(gene) %>%
  unique()

# PDF of all cells

mrk_plot_pos <- de_res %>%
  dplyr::filter(gene %in% pos_genes) %>%
  tidyr::complete(name, gene, fill = list(logFC = 0)) %>%
  ggplot(aes(x = name, y = factor(gene,
                                  levels = pos_genes), fill = logFC)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
        axis.text.y = element_text(size = 7)) +
  ylab("") +
  scale_fill_gradient2()

pdf("./results/niche_mapping/niche_up_genes.pdf", height = 8, width = 3.5)
plot(mrk_plot_pos)
dev.off()

# Get top 10 genes per cell
neg_genes <- de_res %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  dplyr::filter(FDR  <= 0.25) %>%
  arrange(name, logFC) %>%
  group_by(name) %>%
  dplyr::slice(1:10) %>%
  pull(gene) %>%
  unique()

# PDF of all cells

mrk_plot_neg <- de_res %>%
  dplyr::filter(gene %in% neg_genes) %>%
  tidyr::complete(name, gene, fill = list(logFC = 0)) %>%
  ggplot(aes(x = name, y = factor(gene,
                                  levels = neg_genes), fill = logFC)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
        axis.text.y = element_text(size = 7)) +
  ylab("") +
  scale_fill_gradient2()

pdf("./results/niche_mapping/niche_down_genes.pdf", height = 8, width = 3.5)
plot(mrk_plot_neg)
dev.off()

# RUN correlation of niche expression

all_useful_genes <- de_res %>%
  arrange(name, PValue) %>%
  group_by(name) %>%
  dplyr::slice(1:100) %>%
  pull(gene) %>%
  unique()

signature_cor <- de_res %>%
  dplyr::filter(gene %in% all_useful_genes) %>%
  tidyr::complete(name, gene, fill = list(logFC = 0)) %>%
  dplyr::select(name,gene, logFC) %>%
  pivot_wider(names_from = name, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column("niche_a") %>%
  pivot_longer(-niche_a, names_to = "niche_b", values_to = "corr")

signature_cor_plt <- ggplot(signature_cor, aes(x = niche_a, y = niche_b, fill = corr)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_equal() +
  scale_fill_gradient2(lim = c(-1,1))

pdf("./results/niche_mapping/correlation_niche_footprints.pdf", height = 3, width = 3)
plot(signature_cor_plt)
dev.off()

# Run PROGENy

gls <- dplyr::select(de_res, 
                     name, gene, 
                     PValue, logFC) %>%
  dplyr::mutate(pval = -log10(PValue)) %>%
  mutate(pval = ifelse(is.infinite(pval), NA, pval))

max_pval <- max(gls$pval,na.rm = T)

gls <- gls %>%
  dplyr::mutate(pval = ifelse(is.na(pval), max_pval, pval)) %>%
  dplyr::mutate(rank_ord = pval * logFC)

progeny_res <- gls %>% 
  dplyr::select(-c("pval", "PValue", "logFC")) %>%
  pivot_wider(names_from = name,
              values_from = rank_ord,
              values_fill = 0) %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  progeny(.,scale = F,top = 500,perm = 1000)


progeny_order <- hclust(dist(t(progeny_res)))
progeny_order <- progeny_order$labels[progeny_order$order]

progeny_plt <- progeny_res %>%
  as.data.frame() %>%
  rownames_to_column("niche") %>%
  pivot_longer(-niche, names_to = "pathway", values_to = "progeny_scores") %>%
  mutate(pathway = factor(pathway, levels = progeny_order)) %>%
  ggplot(aes(x = pathway, y = niche, fill = progeny_scores)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_equal()


pdf("./results/niche_mapping/niche_progeny.pdf", height = 4, width = 5)
plot(progeny_plt)
dev.off()

# Show number of genes per niche

de_res %>%  
  dplyr::filter(logFC > 0.5, PValue < 0.05) %>%
  group_by(name) %>%
  summarise(n_genes_up = n())


