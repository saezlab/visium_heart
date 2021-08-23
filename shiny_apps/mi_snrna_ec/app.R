#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(SingleCellExperiment)
library(iSEE)

stopifnot(packageVersion("iSEE") >= "2.4.0")

ec_data <- readRDS("testcells.rds")

# Left panel: Reduced dimensions (UMAP)
panel_reddim <- ReducedDimensionPlot(PanelWidth=6L, Type = "UMAP")
panel_reddim[["ColorBy"]] = "Column data"
panel_reddim[["ColorByColumnData"]] <- "opt_state"
panel_reddim[["Downsample"]] <- TRUE

# Right panel: Violin plots
panel_feats <- FeatureAssayPlot(PanelWidth=6L)
panel_feats[["FacetByColumn"]] <- "opt_state"
panel_feats[["Downsample"]] <- TRUE

#tour <- read.csv(sep = "\t", 
#                 file = "ec_tour.txt", 
#                 comment.char = "")

app <- iSEE(ec_data, 
            appTitle = "Kuppe, Ramirez Flores, Li, et al. Test Cells",
            initial=list(
                panel_reddim,
                panel_feats))#,
            #tour = tour)

shiny::runApp(app)
