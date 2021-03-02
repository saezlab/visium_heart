library(stringr)
library(biomaRt)

mouse_genes <- read.csv("./markers/mousegenes.csv",
                        header = F, 
                        stringsAsFactors = F)[,1]

de_genes <- readRDS("./visium_results_manuscript/ct_data/cardio_coembd/cardio_dea.rds")

de_genes_f <- readRDS("./visium_results_manuscript/ct_data/fibro_coembd/fibro_dea.rds")


# First load everything
mouse_genes <- unlist(map(mouse_genes, function(x) {
  strsplit(split = "/",x)
}))

#' Basic function to convert human to mouse gene names
#' @param x: a vector of genes to be transformed 
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"),
                   martL = human, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  return(genesV2)
}

# Convert everything
converted_genes <- convertHumanGeneList(mouse_genes)

# Check for genes not captured
not_captured_genes <- mouse_genes[!mouse_genes %in% converted_genes$MGI.symbol]
not_captured_genes <- tolower(not_captured_genes)
not_captured_genes <- str_to_title(not_captured_genes)

converted_genes_bis <- convertHumanGeneList(not_captured_genes)

converted_genes = rbind(converted_genes, converted_genes_bis)


# Check

de_genes$RNA %>%
  dplyr::filter(p_val_adj <= 0.01,
                gene %in% converted_genes$HGNC.symbol) %>%
  arrange(cluster, -avg_logFC) %>%
  write.table("./visium_results_manuscript/ct_data/cardio_coembd/cardio_gene_intersection.txt",
              col.names = T,
              sep = "\t",
              row.names = F,
              quote = F)

de_genes$dorothea %>%
  dplyr::filter(p_val_adj <= 0.01,
                gene %in% converted_genes$HGNC.symbol) %>%
  arrange(cluster, -avg_logFC) %>%
  write.table("./visium_results_manuscript/ct_data/cardio_coembd/cardio_TFdorothea_intersection.txt",
              col.names = T,
              sep = "\t",
              row.names = F,
              quote = F)

de_genes_f$RNA %>%
  dplyr::filter(p_val_adj <= 0.01,
                gene %in% converted_genes$HGNC.symbol) %>%
  arrange(cluster, -avg_log2FC) %>%
  write.table("./visium_results_manuscript/ct_data/fibro_coembd/fibro_gene_intersection.txt",
              col.names = T,
              sep = "\t",
              row.names = F,
              quote = F)

de_genes_f$dorothea %>%
  dplyr::filter(p_val_adj <= 0.01,
                gene %in% converted_genes$HGNC.symbol) %>%
  arrange(cluster, -avg_log2FC) %>%
  write.table("./visium_results_manuscript/ct_data/fibro_coembd/fibro_TFdorothea_intersection.txt",
              col.names = T,
              sep = "\t",
              row.names = F,
              quote = F)

de_genes_f$RNA[de_genes_f$RNA$gene %in% converted_genes$HGNC.symbol, ]

de_genes$RNA[de_genes$RNA$gene %in% converted_genes$HGNC.symbol, ] %>%
  dplyr::fil


library(progeny)

gene_scores = getModel(top =100)

test = gene_scores %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  dplyr::filter(value != 0)


test = tibble(EGFR = c(1,1,1,0,0,0),
       MAPK = c(1,1,1,1,1,1),
       EGFR_pert = c(1,1,1,1,1,1),
       MAPK_pert = c(1,1,1,1,1,1))


lm(EGFR_pert ~ 0 + EGFR + MAPK, 
   data = test)






