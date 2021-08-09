#' Creates a heatmap of markers from any
#' Seurat object and clusters
#' @param SeuratObject = seurat object
#' @param GSC = marker data frame with two columns gene & cell_type
#' @param assay = "RNA", assay to plot 
#' @param identity_label = "cell_type", used to separate cell clusters
#' @param show_hr = TRUE
#' 
library(ComplexHeatmap)

domarker_hmap <- function(SeuratObject, 
                          GSC, 
                          assay = "RNA", 
                          identity_label = "cell_type",
                          show_hr=TRUE) {
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 35, c = 100)[1:n]
  }
  
  mat <- SeuratObject@assays[[assay]]@scale.data
  
  GSC = GSC[GSC[, 1] %in% rownames(mat), ]
  
  genes <- GSC[, 1]  
  genes.cols <- GSC[, 2]
  
  cl <- as.character(SeuratObject@meta.data[, identity_label])
  
  # Reorder
  ord <- order((cl), decreasing = FALSE)
  mat <- mat[,ord]
  cl <- cl[ord]
  
  cl.cols <- setNames(gg_color_hue(length(unique(cl))),unique(as.character(sort((cl)))))
  
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
                                                        title_gp=gpar(fontsize=10),
                                                        # at=c(-2.5,-2,-1,0,1,2,2.5),
                                                        # labels = c("","-2","-1","0","1","2",""),
                                                        labels_gp = gpar(fontsize = 10)))
  
  hr <- rowAnnotation(df=data.frame(markers=genes.cols), 
                      show_annotation_name = FALSE,
                      show_legend = TRUE)
  
  f1 <-  circlize::colorRamp2(c(-5, 0, +5), c("purple", "black", "yellow"))
  hp <- Heatmap(mat2, cluster_rows = FALSE, cluster_columns = FALSE,col = f1,
                name="Expression",
                top_annotation = hc, bottom_annotation = hc,
                split=factor(genes.cols, levels=unique(genes.cols)),row_title_rot = 0,row_gap = unit(1.5, "mm"), 
                column_split = factor(cl, levels=unique(cl)), column_title_rot=90, column_gap = unit(0.7, "mm"),
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