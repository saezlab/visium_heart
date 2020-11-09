library(Matrix)
library(Seurat)
library(optparse)
library(RANN)
library(data.table)
library(doParallel)
registerDoParallel(cores=20)

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
    parser <- add_option(parser, c("-s", "--SAMPLE"), type="character", default="CK166",
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


#coembed <-readRDS(file=paste0("save/transfered_merged_data_", atac_name , ".Rds"))
#rna_a <- readRDS(file=paste0("save/", rna_name, ".Rds"))
#atac <- readRDS(file=paste0("save/", atac_name, ".Rds"))
atac_dir <- paste0("../ATAC_Single_Sample_V2/data/", atac_name)
atac <- readRDS(file=file.path(atac_dir,paste0(atac_name, "_unionPeaks.Rds")))

DefaultAssay(atac) <- "peaks"


a <- colnames(atac)
b <- colnames(atac)

rst <- FindNN(atac, a, b, reduction="lsi", k=60, dims=2:30)
nns <- rst@tools$integrated@neighbors

rm(rst)
gc()


peak_nearest_cell_list = list()
for(i in seq_along(colnames(atac))){
  vec <- b[nns$nnbb$nn.idx[i, 1:51]]
  peak_nearest_cell_list[[i]] <- vec 
}



fname = paste0("save/", atac_name, "_atac_aggre.csv")
write.table(transpose(peak_nearest_cell_list), sep = "\t", 
            file = fname, row.names = T,col.names = F,
            quote=FALSE)


cell_num <- length(peak_nearest_cell_list)
#atac <- readRDS(file=paste0("save/", atac_name, ".Rds"))
mtx_atac <- GetAssayData(atac, slot="counts", assay="peaks")
aggre_atac_mtx <- foreach(i =1:cell_num,  .combine = cbind)%dopar%{
                one <- peak_nearest_cell_list[[i]]
                cell_name <- one[1]
                rowSums(mtx_atac[, one[1:51]])
}
colnames(aggre_atac_mtx) <- colnames(atac)
print(paste(atac_name,"end", date()))
saveRDS(aggre_atac_mtx, file=sprintf("save/%s_aggred_atac.Rds", atac_name))

