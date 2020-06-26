library(Seurat)
library(tidyverse)
library(cowplot)
library(readbitmap)
library(grid)
library(rjson)


sample_names <- c("AKK006 - healthy control", "AKK004 - old MI, fibrotic areas", 
                  "AKK003 - acute MI", "AKK003 - borderzone", # "AKK002 - acute MI", 
                  "AKK002 - borderzone", "AKK002 - healthy part", 
                  "AKK001 - fibrotic, late stage heart failure")

matrix_paths <- c("157771/157771/outs/filtered_feature_bc_matrix.h5",
                  "157772/157772/outs/filtered_feature_bc_matrix.h5",
                  "157775/157775/outs/filtered_feature_bc_matrix.h5",
                  "157777/157777/outs/filtered_feature_bc_matrix.h5",
                  #"157779/157779/outs/filtered_feature_bc_matrix.h5",
                  "157781/157781/outs/filtered_feature_bc_matrix.h5",
                  "157782/157782/outs/filtered_feature_bc_matrix.h5",
                  "157785/157785/outs/filtered_feature_bc_matrix.h5")


integrated_data <- c("../ATAC_RNA_Integration/integrated/CK166/CK166_integrated_unionPeaks.Rds",
                     "../ATAC_RNA_Integration/integrated/CK170/CK170_integrated_unionPeaks.Rds",
                     "../ATAC_RNA_Integration/integrated/CK174/CK174_integrated_unionPeaks.Rds",
                     "../ATAC_RNA_Integration/integrated/CK171/CK171_integrated_unionPeaks.Rds",
                     # "../ATAC_RNA_Integration/integrated/CK169/CK169.Rds",
                     "../ATAC_RNA_Integration/integrated/CK168/CK168_integrated_unionPeaks.Rds",
                     "../ATAC_RNA_Integration/integrated/CK173/CK173_integrated_unionPeaks.Rds",
                     "../ATAC_RNA_Integration/integrated/CK167/CK167_integrated_unionPeaks.Rds")

rnaseq_data <- c("../ATAC_RNA_Integration/RNA/CK158.filtered.annotated.Rds",
                 "../ATAC_RNA_Integration/RNA/CK162.filtered.annotated.Rds",
                 "../ATAC_RNA_Integration/RNA/CK165.filtered.annotated.Rds",
                 "../ATAC_RNA_Integration/RNA/CK163.filtered.annotated.Rds",
                 #"../ATAC_RNA_Integration/RNA/CK161.Rds",
                 "../ATAC_RNA_Integration/RNA/CK160.filtered.annotated.Rds",
                 "../ATAC_RNA_Integration/RNA/CK164.filtered.annotated.Rds",
                 "../ATAC_RNA_Integration/RNA/CK159.filtered.annotated.Rds"
)


seq_along(sample_names) %>% walk(function(i){
  
  spatial <- CreateSeuratObject(Read10X_h5(matrix_paths[i]), project = sample_names[i]) %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
  
  integrated <- readRDS(integrated_data[i])
  
  # we need the rna to get the variable features that were used 
  # for the integration of rnaseq and atacseq
  rna <- readRDS(rnaseq_data[i])

  anchors <- FindTransferAnchors(reference = integrated, reference.assay = "RNA", 
                                 features = rna@assays$RNA@var.features, query = spatial, 
                                 reduction = "cca")
  
  predictions.assay <- TransferData(anchorset = anchors, refdata = integrated$celltype, 
                                    weight.reduction = spatial[["pca"]])
  
  spatial <- AddMetaData(spatial, metadata = predictions.assay)
  
  write_rds(spatial, path = paste0("transfer_results/integrated/", str_split(matrix_paths[i],"/")[[1]][1], ".rds"))
})

metas <- list.files("transfer_results/integrated/",full.names = T) %>% map(function(f){
  d <- read_rds(f)
  d@meta.data
})
names(metas) <- str_split(list.files("transfer_results/integrated/"),"\\.") %>% map(~.x[1])
write_rds(metas, "transfer_results/integrated_meta.rds")

geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- viewport(x=data$x, y=data$y)
      g <- editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))



image_paths <- c("157771/157771/outs/spatial/tissue_hires_image.png",
                 "157772/157772/outs/spatial/tissue_hires_image.png",
                 "157775/157775/outs/spatial/tissue_hires_image.png",
                 "157777/157777/outs/spatial/tissue_hires_image.png",
#                 "157779/157779/outs/spatial/tissue_hires_image.png",
                 "157781/157781/outs/spatial/tissue_hires_image.png",
                 "157782/157782/outs/spatial/tissue_hires_image.png",
                 "157785/157785/outs/spatial/tissue_hires_image.png")


scalefactor_paths <-c("157771/157771/outs/spatial/scalefactors_json.json",
                      "157772/157772/outs/spatial/scalefactors_json.json",
                      "157775/157775/outs/spatial/scalefactors_json.json",
                      "157777/157777/outs/spatial/scalefactors_json.json",
#                      "157779/157779/outs/spatial/scalefactors_json.json",
                      "157781/157781/outs/spatial/scalefactors_json.json",
                      "157782/157782/outs/spatial/scalefactors_json.json",
                      "157785/157785/outs/spatial/scalefactors_json.json")

tissue_paths <- c("157771/157771/outs/spatial/tissue_positions_list.txt",
                  "157772/157772/outs/spatial/tissue_positions_list.txt",
                  "157775/157775/outs/spatial/tissue_positions_list.txt",
                  "157777/157777/outs/spatial/tissue_positions_list.txt",
#                  "157779/157779/outs/spatial/tissue_positions_list.txt",
                  "157781/157781/outs/spatial/tissue_positions_list.txt",
                  "157782/157782/outs/spatial/tissue_positions_list.txt",
                  "157785/157785/outs/spatial/tissue_positions_list.txt")


images_cl <- list()

for (i in 1:length(sample_names)) {
  images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
  height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
  width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)

grobs <- list()
for (i in 1:length(sample_names)) {
  grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width

scales <- list()

for (i in 1:length(sample_names)) {
  scales[[i]] <- fromJSON(file = scalefactor_paths[i])
}

bcs <- list()

for (i in 1:length(sample_names)) {
  bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
  bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_hires_scalef    # scale tissue coordinates for lowres image
  bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_hires_scalef
  bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
  bcs[[i]]$height <- height$height[i]
  bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names


bcs_merged <- bcs %>% map2(metas, function(x,y){
  merge(x, y%>%rownames_to_column("barcode"))
})

plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merged[[i]] %>% filter(prediction.score.max > 0.5) %>%
    ggplot(aes_string(x="imagecol",y="imagerow", fill= "predicted.id")) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 0.8, stroke = 0.3)+
    coord_cartesian(expand=FALSE)+    xlim(c(0,max(bcs_merged[[i]] %>%
                   select(width))))+
    ylim(c(max(bcs_merged[[i]] %>%
                 select(height)),0))+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

pdf("transfer_results/labels_integrated.pdf", width=15, height=10)
plot_grid(plotlist = plots, align = "v")
dev.off()

