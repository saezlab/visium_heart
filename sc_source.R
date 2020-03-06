# Find inflection points
get_npcs = function(seurat_object, show_plot = T){
  library(ggplot2)
  
  std_data = Stdev(object = seurat_object, reduction = "pca")
  ndims = length(x = std_data)
  elbow_data = data.frame(dims = 1:ndims, stdev = std_data[1:ndims])
  
  reference = elbow_data[,2]
  difference_vect = c(elbow_data[2:nrow(elbow_data),2],0)
  difference = which((reference - difference_vect) < 0.05)
  
  difference_consec = c(difference[2:length(difference)],0) - difference
  names(difference_consec) = difference
    
  npcs = as.numeric(names(which(difference_consec ==1)[1]))
  
  if(show_plot){
    
    plt = ggplot(elbow_data, aes(x = dims,
                                 y = stdev)) +
          geom_point() + geom_vline(xintercept=npcs) +
          theme_minimal()
    
    plot(plt)
  }
  
  return(npcs)
  
}

# Enhanced DoHeatmap function
DoHeatmap2 <- function(SeuratObject, GSC, assay="RNA", res=0.5, show_hr=TRUE) {
  library(ComplexHeatmap)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  GSC = GSC[GSC[,1] %in% rownames(mat),]
  
  genes = GSC[,1]  
  genes.cols = GSC[,2]
  
  if(is.null(res)) {
    cl <- as.character(SeuratObject@meta.data[,"seurat_clusters"])
  } else {
    cl <- as.character(SeuratObject@meta.data[,paste0(assay,"_snn_res.",res)]) 
  }
  
  # Reorder
  ord <- order(as.numeric(cl), decreasing = FALSE)
  mat <- mat[,ord]
  cl <- cl[ord]
  
  cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort(as.numeric(cl)))))
  
  common_genes <- intersect(genes,rownames(mat))
  diff_genes <- setdiff(genes, rownames(mat))
  
  mat2 <- rbind(mat[common_genes,],
                matrix(NA, nrow=length(diff_genes), ncol=ncol(mat),
                       dimnames=list(diff_genes, colnames(mat)))
  )
  mat2 <- mat2[genes,]
  
  hc <- HeatmapAnnotation(df=data.frame("cluster"=cl),col = list("cluster"=cl.cols), show_annotation_name = FALSE,
                          show_legend = FALSE,
                          annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                        grid_width = unit(5, "mm"),
                                                        title_gp=gpar(fontsize=16),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 14)))
  
  hr <- rowAnnotation(df=data.frame(markers=genes.cols), 
                      show_annotation_name = FALSE,
                      show_legend = TRUE)
  
  f1 <-  circlize::colorRamp2(c(-5, 0, +5), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name="Expression",
                top_annotation = hc, bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_gap = unit(0.7, "mm"),
                show_column_names = FALSE, row_names_side = "left",
                show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize=7))
  if(show_hr) {
    hh <- hp + hr 
  } else {
    hh <- hp
  }
  return(hh)
}



## Function to run viper combined with Dorothea Regulons on single cell data.
## Inputs:
## - InputObject: An object of class Seurat or SingleCellExperiment. 
## The normalized expression values will be extracted from the object.
## It can also be a matrix containing normalized expresion values.
## - regulon: Object of class regulon. See viper package.
## - ... : Additional parameters for Viper function.
## Output:
## - A matrix containing the activity of the different TFs provided in the 
## regulon object.
sc_viper = function(InputObject, regulon, ...) {
  
  if (class(InputObject) == "Seurat"){
    expr <- InputObject[["RNA"]]@data
  } else {
    if(class(InputObject) == "SingleCellExperiment"){
      expr <- normcounts(InputObject)
    } else {
      expr <- InputObject
    }
  }
  
  res <- viper(eset = as.matrix(expr), regulon = regulon, nes = TRUE, 
               method = "scale", minsize = 4, eset.filter = FALSE, 
               cores = 1, verbose = FALSE)
  return(res)
}

## Function to group Dorothea regulons. 
## Input: A data frame containing Dorothea regulons, as stored in 
## https://github.com/saezlab/ConservedFootprints/tree/master/data
## Output: Object of class regulon. See viper package.
df2regulon = function(df) {
  regulon = df %>%
    split(.$tf) %>%
    map(function(dat) {
      tf = dat %>% distinct(tf) %>% pull()
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode =targets, likelihood = likelihood)
    })
  return(regulon)
}


#' Function to calculate differentially expressed genes
#' 
#' Inputs:
#' -Seurat object
#' -diseases: a vector len = 2 with the conditions to compare
#' -cluster_ids = named vector with idents to test independently 
#'  
#' Outputs:
#' A list length = length(cluster_ids) with results from the test
#' 
get_deg = function(seurat_obj, diseases, cluster_ids){
  
  cluster_deg = lapply(cluster_ids, function(c){
    
    contrast_vect = paste(c,diseases,sep = "_")
    
    deg = FindMarkers(seurat_obj, 
                      ident.1 = contrast_vect[1], 
                      ident.2 = contrast_vect[2], 
                      verbose = FALSE,
                      logfc.threshold = 0.01
    )
    
    return(deg)
    
  })
  
  return(cluster_deg)
  
}


#' Generates hmaps of marker TFs/genes for publication 
#' @param mat Scaled expression matrix. Features in rows, samples in columns. Must be ordered based on an identity, like cell type
#' @param GSC Data frame of gene markers, genes in first column, cluster on the right
#' @param col_palette color pallete size of numbers of identities in mat
#' @return complex heatmap plot object
#' 
plot_scpaperhmap = function(mat,GSC,col_palette = col_palette){
  require(ComplexHeatmap)
  # Assign colors to classes
  cols_df = data.frame(group_name = unique(cl),
                       colorname = col_palette,
                       stringsAsFactors = F)
  
  cl.cols = col_palette
  cl.cols = left_join(data.frame(group_name=cl,stringsAsFactors = F),
                      cols_df)[,"colorname"]
  
  names(cl.cols) = cl
  
  # Upper 
  hc = HeatmapAnnotation(df=data.frame("cluster"=cl),
                         col = list("cluster"=cl.cols), 
                         show_annotation_name = FALSE,
                         show_legend = TRUE,
                         annotation_legend_param= list(legend_height = unit(8, "cm"),
                                                       grid_width = unit(5, "mm"),
                                                       title_gp=gpar(fontsize=10),
                                                       # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                       # labels = c("","-2","-1","0","1","2",""),
                                                       labels_gp = gpar(fontsize = 12)))
  
  f1 =  circlize::colorRamp2(c(-2,0,+2), c("purple", "black", "yellow"))
  
  hp = Heatmap(mat, 
               cluster_rows = FALSE, 
               cluster_columns = FALSE,col = f1,
               name="Expression",
               top_annotation = hc, bottom_annotation = NULL,
               #split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
               #column_split = factor(cl, levels=unique(cl)), column_title_rot=00, column_gap = unit(0.7, "mm"),
               show_column_names = FALSE, 
               row_names_side = "left",
               show_heatmap_legend = TRUE,
               row_names_gp = gpar(fontsize=12))
  
  
  return(hp)
}