# Copyright (c) [2021] [Ricardo O. Ramirez Flores, Jesus Velez]
# roramirezf@uni-heidelberg.de

#' Run and deploy spatialLIBD app

library("spatialLIBD")
library("markdown")

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
spe <- readRDS("spe_libd.rds")
other_continous <- colnames(colData(spe))
progeny <- other_continous[grepl(pattern = "progeny", other_continous)]
c2l <- other_continous[grepl(pattern = "c2l", other_continous)]


## Deploy the website
spatialLIBD::run_app(
  spe,
  sce_layer = NULL,
  modeling_results = NULL,
  sig_genes = NULL,
  title = "Specimen AKK001_157785",
  spe_discrete_vars = c("label", 
                        "ManualAnnotation",
                        "opt_clust",
                        "opt_clust_integrated"),
  spe_continuous_vars = c(
    "sum_umi",
    "sum_gene",
    "expr_chrM",
    "expr_chrM_ratio",
    "sum",
    "detected",
    "subsets_mito_sum",
    "subsets_mito_detected",
    "subsets_mito_percent",
    "total",
    "sizeFactor",
    progeny,
    c2l
  ),
  default_cluster = "label",
  docs_path = "www"
)
