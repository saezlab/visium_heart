# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate markers of cells
#' using edgeR and pseudobulk profiles of all samples
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)


pseudobulk_data <- readRDS("./processed_snrnaseq/integration/psxsmpl_integrated_rnasamples_ann.rds")[[1]]
patient_meta <- readRDS("./markers/snrna_patient_anns_revisions.rds")

counts <- assay(pseudobulk_data[["gex"]])
pb_meta <- colData(pseudobulk_data[["gex"]])

pb_meta <- pb_meta %>% 
  as.data.frame() %>% 
  left_join(patient_meta, 
            by = c("orig.ident" = "sample_id"))

patient_cells <- pb_meta %>%
  group_by(patient_id, cell_type) %>%
  summarize(ncells = sum(ncells))


# This summarizes the info into patients
pseudobulk_data <- sumCountsAcrossCells(counts, 
                                DataFrame(pb_meta[, c("patient_id", "cell_type")]))

rm(pb_meta)

counts <- assay(pseudobulk_data)

meta_data <- colData(pseudobulk_data)

meta_data$max_counts <- colMaxs(counts)

meta_data <- meta_data %>%
  as.data.frame() %>%
  dplyr::select(-ncells) %>%
  left_join(patient_cells)

# We will exclude all profiles that are coming from less than 10 cells and a min of 100 max counts
meta_data$keep <- ifelse(meta_data$ncells < 50 &
                         meta_data$max_counts < 1000, 
                         FALSE, TRUE)

cols <- c("TRUE" = "grey",
          "FALSE" = "black")

# First check the number of cells per patient
patient_summary_plt <- ggplot(meta_data, 
       aes(y = patient_id, x = cell_type, fill = keep)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background= element_rect(fill="black", colour="black")) +
  scale_fill_manual(values = cols,
                     na.value = 'black')

t_ix <- which(meta_data$keep == TRUE)

meta_data <- meta_data[t_ix, ]
counts <- counts[, t_ix]

meta_data$col_id <- paste0(meta_data$patient_id,
                            "_",
                            meta_data$cell_type)

colnames(counts) <- meta_data$col_id 

# From here we loop through CTs

cts <- set_names(meta_data$cell_type %>% unique)

de_res <- map(cts, function(ct) {
  print(ct)
  ct_meta_data <- meta_data %>%
    mutate(test_column = ifelse(cell_type == ct, ct, "rest"))
  
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
  write_csv(file = "./results/cell_markers/edgeR_cellmrkrs.csv")

# Get top 10 genes per cell

genes <- de_res %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  arrange(name, FDR, - logFC) %>%
    group_by(name) %>%
    dplyr::slice(1:10) %>%
  pull(gene)


# PDF of all cells

mrk_plot <- de_res %>%
  dplyr::filter(gene %in% genes) %>%
  tidyr::complete(name, gene, fill = list(FDR = 1)) %>%
  ggplot(aes(x = name, y = factor(gene,
                                  levels = genes), fill = -log10(FDR))) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
        axis.text.y = element_text(size = 7)) +
  ylab("")

pdf("./results/cell_markers/edgeR_cellmrkrs_snRNA.pdf", height = 16, width = 4)  
plot(mrk_plot)
dev.off()

# This is the object of cell_type markers

marker_list <- de_res %>%
  dplyr::filter(logFC > 0, FDR < 0.15) %>%
  arrange(name, FDR, - logFC) %>%
  dplyr::filter(!grepl("^AC[0-9]", gene)) %>%
  dplyr::filter(!grepl("^AL[0-9]", gene)) %>%
  dplyr::filter(!grepl("^LINC[0-9]", gene)) %>%
  arrange(name, FDR, -logFC) %>%
  group_by(name) %>%
  dplyr::slice(1:200) %>%
  dplyr::select(name, gene) %>%
  nest() %>%
  mutate(data = map(data, ~ .x[[1]] %>%
                      sort())) %>%
  deframe()

saveRDS(marker_list, "./markers/pb_ct_marker_list.rds")
