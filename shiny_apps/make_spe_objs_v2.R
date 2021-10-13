# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we create SpatialLIBD ready objects

## "Visium human DLPFC workflow"
library("dplyr")
library("SpatialExperiment")
library("scater")
library("scran")
library("igraph")
library("BiocFileCache")
library("rtracklayer")
library("lobstr")
library("spatialLIBD")
library("Seurat")
library("purrr")

# Specify working directories ---------------------------------------------
visium_data_dir <- "./visium_data/"
seurat_data_dir <- "./processed_visium/objects/"

param_df <- tibble(raw = paste0(visium_data_dir,
                                list.dirs(visium_data_dir,
                                          recursive = F,
                                          full.names = F),
                                "/outs"),
                   processed = paste0(seurat_data_dir,
                                      list.dirs(visium_data_dir,
                                                recursive = F,
                                                full.names = F), 
                                      ".rds"),
                   sample_id = list.dirs(visium_data_dir,
                                         recursive = F,
                                         full.names = F),
) 

raw = "./visium_data/AKK001_157785/outs"
processed = "./processed_visium/objects/AKK001_157785.rds"
sample_id = "AKK001_157785"

# This creates a spe from scratch and copies useful info
# You only have raw data here
create_spe <- function(raw, processed, sample_id) {
  
  print(sample_id)
  spe_out_dir <- paste0("./shiny_apps/visium/", sample_id)
  system(paste0("mkdir ", spe_out_dir))
  
  spe_out <- paste0(spe_out_dir, "/", "spe_libd", ".rds")
  
  # reading spe data
  spe <- SpatialExperiment::read10xVisium(raw, 
                                          sample_id = sample_id,
                                          images = c("lowres", "hires"),
                                          type = "sparse")
  
  # Transform to compatible data for SpatialLIBD
  spe@assays@data$counts <- as(spe@assays@data$counts, "dgTMatrix")
  
  # Load seurat object -----------------------------
  seurat <- readRDS(processed)
  
  # subset to keep only spots over tissue
  spe <- spe[, spatialData(spe)$in_tissue == 1]
  
  # identify mitochondrial genes
  rowData(spe)$gene_name <- rowData(spe)$symbol
  is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
  
  # calculate per-spot QC metrics and store in colData
  spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
  
  # Subset spots from seurat object ----------------
  spe <- spe[, colnames(seurat)]
 
  # quick clustering for pool-based size factors
  set.seed(123)
  qclus <- quickCluster(spe)
  
  # calculate size factors and store in object
  spe <- computeSumFactors(spe, cluster = qclus)
  
  # calculate logcounts (log-transformed normalized counts) and store in object
  spe <- logNormCounts(spe)
  
  # remove mitochondrial genes
  spe <- spe[!is_mito, ]
  
  # fit mean-variance relationship
  dec <- modelGeneVar(spe)
  
  # select top HVGs
  top_hvgs <- getTopHVGs(dec, prop = 0.1)
  
  # compute PCA
  set.seed(123)
  spe <- runPCA(spe, subset_row = top_hvgs)
  
  # compute UMAP on top 50 PCs
  set.seed(123)
  spe <- runUMAP(spe, dimred = "PCA")
  
  # update column names for easier plotting
  colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
  
  # graph-based clustering
  set.seed(123)
  k <- 10
  g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
  g_walk <- igraph::cluster_walktrap(g)
  clus <- g_walk$membership
  
  # store cluster labels in column 'label' in colData
  colLabels(spe) <- factor(clus)
  
  # set gene names as row names for easier plotting
  rownames(spe) <- rowData(spe)$gene_name
  
  # test for marker genes
  markers <- findMarkers(spe, test = "binom", direction = "up")
  
  ## Find the interesting markers for each cluster
  interesting <- sapply(markers, function(x) x$Top <= 5)
  colnames(interesting) <- paste0("gene_interest_", seq_len(length(markers)))
  rowData(spe) <- cbind(rowData(spe), interesting)
  
  
  spe$key <- paste0(spe$sample_id, "_", colnames(spe))
  spe$sum_umi <- colSums(counts(spe))
  spe$sum_gene <- colSums(counts(spe) > 0)
  
  ## Download the Gencode v32 GTF file and cache it
  bfc <- BiocFileCache::BiocFileCache()
  gtf_cache <- BiocFileCache::bfcrpath(
    bfc,
    paste0(
      "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
      "release_32/gencode.v32.annotation.gtf.gz"
    )
  )
  
  ## Import into R (takes ~1 min)
  gtf <- rtracklayer::import(gtf_cache)
  
  ## Subset to genes only
  gtf <- gtf[gtf$type == "gene"]
  
  ## Remove the .x part of the gene IDs
  gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)
  
  ## Set the names to be the gene IDs
  names(gtf) <- gtf$gene_id
  
  ## Match the genes
  match_genes <- match(rowData(spe)$gene_name, gtf$gene_name)

  ## Drop the few genes for which we don't have information
  spe <- spe[!is.na(match_genes), ]
  match_genes <- match_genes[!is.na(match_genes)]
  
  ## Keep only some columns from the gtf
  mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]
  
  ## Save the "interest"ing columns from our original spe object
  interesting <- rowData(spe)[, grepl("interest", colnames(rowData(spe)))]
  
  ## Add the gene info to our SPE object
  rowRanges(spe) <- gtf[match_genes]
  
  ## Add back the "interest" coolumns
  rowData(spe) <- cbind(rowData(spe), interesting)
  
  ## Add information used by spatialLIBD
  rowData(spe)$gene_search <- paste0(
    rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
  )
  
  ## Compute chrM expression and chrM expression ratio
  is_mito <- which(seqnames(spe) == "chrM")
  spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
  spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi
  
  ## Add a variable for saving the manual annotations
  spe$ManualAnnotation <- "NA"
  
  ## Remove genes with no data
  no_expr <- which(rowSums(counts(spe)) == 0)
  
  spe <- spe[-no_expr, , drop = FALSE]
  
  if (any(spe$sum_umi == 0)) {
    spots_no_counts <- which(spe$sum_umi == 0)
    ## Number of spots with no counts
    print(length(spots_no_counts))
    ## Percent of spots with no counts
    print(length(spots_no_counts) / ncol(spe) * 100)
    spe <- spe[, -spots_no_counts, drop = FALSE]
  }
  
  spatialData(spe) <- cbind(spatialData(spe), barcode_id = rownames(spatialData(spe)))
  
  # Just do some basic adjustments on the spatial coordinates
  
  spatialData(spe) <- cbind(spatialData(spe), spatialCoords(spe))
  spatialCoords(spe) <- as.matrix(DataFrame(x = spatialCoords(spe)[, "pxl_row_in_fullres"],
                                            y = spatialCoords(spe)[, "pxl_col_in_fullres"],
                                            row.names = rownames(spatialCoords(spe))))
  
  # Finally, just copy the data from assays and spots that were kept
  seurat <- seurat[, colnames(spe)]
  
  niche_clustering <- seurat@meta.data[, c("opt_clust", "opt_clust_integrated")] %>% 
    mutate_all(as.character) %>% 
    mutate_all(as.numeric) %>% 
    mutate_all(~ .x + 1)
  
  extra_variables <- get_assay_data(seurat)
  
  colData(spe) <- cbind(
    colData(spe),
    niche_clustering,
    extra_variables
  )
  
  # Save the object
  
  saveRDS(spe, file = spe_out)
  
}

# This function copies the info from the seurat object
# This excludes gene expression data and dorothea
get_assay_data <- function(seurat){
  
  seurat@assays %>% #Exclude SCT and Spatial Assays
    purrr::list_modify(
      Spatial = purrr::zap(),
      SCT = purrr::zap(),
      dorothea = purrr::zap()
    ) %>% 
    purrr::discard(purrr::is_empty) %>%
    purrr::map2(
      .y = names(.),
      .f = ~{
        .x@data %>% 
          `rownames<-`(paste0(.y, "_", rownames(.)))
      }
    ) %>%
    purrr::reduce(base::rbind) %>% 
    `colnames<-`(
      stringr::str_remove(
        string = colnames(.),
        pattern = ".+_"
      )
    ) %>% 
    t() %>% 
    base::identity()
  
}

pwalk(param_df)