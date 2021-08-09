library(scater)
library(HDF5Array)
library(shiny)
library(tidyverse)
library(iSEE)

ec_data <- loadHDF5SummarizedExperiment("./ct_data_sce/integrated_data_wstates/")

# Left panel: Reduced dimensions (UMAP)
panel_reddim <- ReducedDimensionPlot(PanelWidth=6L, Type = "UMAP")
panel_reddim[["ColorBy"]] = "Column data"
panel_reddim[["ColorByColumnData"]] <- "cell_type"
panel_reddim[["Downsample"]] <- TRUE

# Right panel: Violin plots
panel_feats <- FeatureAssayPlot(PanelWidth=6L)
panel_feats[["FacetByColumn"]] <- "cell_type"
panel_feats[["Downsample"]] <- TRUE

tour <- read.csv(sep = "\t", file = "./tours/ec_tour.txt", comment.char = "")

app <- iSEE(ec_data, initial=list(
  panel_reddim,
  panel_feats),
  tour = tour)

shiny::runApp(app)

tour <- defaultTour()




