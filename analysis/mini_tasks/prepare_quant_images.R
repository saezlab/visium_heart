# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Create images for nuclei quantification
#' 
#' Code is adapted from Leo Collado's spatialLIBD
#' https://github.com/LieberInstitute/spatialLIBD/blob/master/R/vis_gene_p.R

library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(ggplot2)

# This is the ggplot2 handler coded by Leo
geom_spatial <- function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  ## To avoid a NOTE on R CMD check
  ggname <- function(prefix, grob) {
    grob$name <- grid::grobName(grob, prefix)
    grob
  }
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x = data$x, y = data$y)
      g <- grid::editGrob(data$grob[[1]], vp = vp)
      ggname("geom_spatial", g)
    },
    required_aes = c("grob", "x", "y")
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


# Function to create samples

plot_rois <- function(sample_path, sample_name) {
  
  spe <- read10xVisium(sample_path,
                       images = "hires")
  
  sampleid <- "sample1"
  ## Prepare the data for the plotting function
  d <- as.data.frame(SpatialExperiment::spatialData(spe, cd_bind = TRUE))
  ## Some variables
  pxl_row_in_fullres <- pxl_col_in_fullres <- key <- COUNT <- NULL
  img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid)
  
  
  # Basis of the plot
  p <-
    ggplot(
      d,
      aes(
        x = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid),
        y = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid),
        fill = NULL,
        color = NULL,
        key = key
      )
    )
  
  # Getting the image
  grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), 
                           height = grid::unit(1, "npc"))
  
  # Generating the empty plot
  p <- p + geom_spatial(
    data = tibble::tibble(grob = list(grob)),
    aes(grob = grob),
    x = 0.5,
    y = 0.5
  ) + geom_point(
    shape = 21,
    size = 12, #3.5
    stroke = 1.5
  ) +
    coord_cartesian(expand = FALSE) +
    xlim(0, nrow(img)) +
    ylim(ncol(img), 0) +
    xlab("") + ylab("") +
    labs(fill = NULL, color = NULL) +
    theme_set(theme_bw(base_size = 20)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = c(0.95, 0.10)
    )
  
  #pdf(width = 16, height = 16, 
  #    file = "./visium_results_manuscript//spot_plots/test.pdf")
  
  #print(p)
  
  #dev.off()
  
  jpeg(file = paste0("./visium_results_manuscript/spot_plots/", 
                     sample_name,
                     ".jpeg"), quality = 100,
       width = 4000, height = 4000, pointsize = 50)
  
  print(p)
  
  dev.off()
  
  
}

# Main
param_df <- tibble(sample_path = paste0(list.dirs("./visium_data", recursive = F), "/"),
                   sample_name = list.files("./visium_data/"))

pmap(param_df, plot_rois)
