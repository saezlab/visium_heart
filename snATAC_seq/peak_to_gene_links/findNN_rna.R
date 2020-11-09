library(Matrix)
library(Seurat)
library(optparse)
library(RANN)
library(data.table)
library(doParallel)
registerDoParallel(cores=1)

FindNN <- function(
  object,
  cells1 = NULL,
  cells2 = NULL,
  grouping.var = NULL,
  dims = 1:10,
  reduction = "cca.l2",
  nn.dims = dims,
  nn.reduction = reduction,
  k = 300, 
  eps = 0, 
  integration.name = 'integrated',
  verbose = TRUE 
) {
  if (xor(x = is.null(x = cells1), y = is.null(x = cells2))) {
    stop("cells1 and cells2 must both be specified")
  }
  if (!is.null(x = cells1) && !is.null(x = cells2) && !is.null(x = grouping.var)) {
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if (is.null(x = cells1) && is.null(x = cells2) && is.null(x = grouping.var)) {
    stop("Please set either cells1/2 or grouping.var")
  }
  if (!is.null(x = grouping.var)) {
    if (nrow(x = unique(x = object[[grouping.var]])) != 2) { 
      stop("Number of groups in grouping.var not equal to 2.")
    }    
    groups <- names(x = sort(x = table(... = object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  dim.data.self <- Embeddings(object = object[[nn.reduction]])[ ,nn.dims]
  dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
  dims.cells1.self <- dim.data.self[cells1, ]
  dims.cells1.opposite <- dim.data.opposite[cells1, ]
  dims.cells2.self <- dim.data.self[cells2, ]
  dims.cells2.opposite <- dim.data.opposite[cells2, ]

  nnaa <- nn2( 
    data = dims.cells1.self,
    k = k + 1, 
    eps = eps
  )
  nnab <- nn2( 
    data = dims.cells2.opposite,
    query = dims.cells1.opposite,
    k = k, 
    eps = eps
  )
  nnbb <- nn2( 
    data = dims.cells2.self,
    k = k + 1, 
    eps = eps
  )
  nnba <- nn2( 
    data = dims.cells1.opposite,
    query = dims.cells2.opposite,
    k = k, 
    eps = eps
  )
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  return(object)
}


AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-s", "--SAMPLE"), type="character", default="CK174",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)
atac_name  = pa$SAMPLE
atac_to_rna = setNames( c("CK158", "CK159", "CK160",  "CK162", "CK163", "CK164","CK165"), 
                    c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174"))
rna_name  = atac_to_rna[atac_name]


atac_dir <- paste0("../ATAC_Single_Sample_V2/data/", atac_name)


rna_dir <- "../scRNA_filtered/"

coembed <-readRDS(file=file.path(atac_dir,paste0(atac_name, "_integrated_unionPeaks.Rds")))
rna_a <- readRDS(file=file.path(rna_dir, sprintf("%s.filtered.annotated.rds", rna_name)))
atac_a <- readRDS(file=file.path(atac_dir, paste0(atac_name, "_unionPeaks.Rds")))

rna_a$celltype = as.character(Idents(rna_a))
rna_a <- subset(rna_a, celltype %in% atac_a$celltype)
saveRDS(rna_a, file=paste0("save/", rna_name, ".Rds"))

a <- colnames(coembed[, coembed$tech=='RNA'])
b <- colnames(coembed[, coembed$tech=='ATAC'])

a <- intersect(a, colnames(rna_a))

rst <- FindNN(coembed, a, b, reduction="pca", k=60)
nns <- rst@tools$integrated@neighbors

rm(rst)
gc()

### RNA for cells reorder
rna_nearest_cell_list = list()
nms <- colnames(atac_a)
for(i in seq_along(colnames(atac_a))){
  #vec <- c(nms[i], a[nns$nnba$nn.idx[i, 1]])
  rna_nm <- a[nns$nnba$nn.idx[i, 1]]
  vec <- c(nms[i],atac_a$celltype[nms[i]], rna_nm, rna_a$celltype[rna_nm])
  rna_nearest_cell_list[[i]] <- vec
}

rm(atac_a)
rm(coembed)
gc()

fname = paste0("save/rna_nearest_",rna_name, ".csv")
write.table(transpose(rna_nearest_cell_list), sep = "\t", 
            file = fname, row.names = T,col.names = F,
            quote=FALSE)

fname = paste0("save/rna_counts_", rna_name, ".txt")
#------write.table(as.matrix(rna_a@assays$RNA@counts), file=fname, quote=FALSE)
write.table(rna_a@assays$RNA@counts, file=fname, quote=FALSE)


counts <- rna_a@assays$RNA@counts
tiny_df = as.data.frame(counts[, 1, drop=F])
#build an empty df
df = tiny_df[,FALSE]

len_cells <- length(rna_nearest_cell_list)
df <- foreach (i = 1:len_cells, .combine=cbind) %do% { 
     counts[, rna_nearest_cell_list[[i]][3]]
     #colnames(df)[i] <- rna_nearest_cell_list[[i]][3]
     #df
}
for(i in 1:len_cells){
    #df[i, ] <- counts[, rna_nearest_cell_list[[i]][3]]
    colnames(df)[i] <- rna_nearest_cell_list[[i]][3]
}
rownames(df) <- rownames(counts)

saveRDS(as.matrix(df), file=sprintf("save/aggred_counterpart_%s.Rds", rna_name))
