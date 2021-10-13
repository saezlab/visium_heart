
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spatial Multiomics Map of Myocardial Infarction

This website allows to visualize gene expression, pathway activities and cell-type annotations of one slide of 10x Visium spatial transcriptomics. 

Complete documentation of this app can be found in the spot-level tab 

For access to objects used in the original publication, please refer to:

## Data transformations 

For visualization purposes, single slides were processed as suggested in the _Visium spatialLIBD workflow_ chapter from  _Orchestrating Spatial Transcriptomics Analyses_ (OSTA) with Bioconductor.

For the same reasons, PCA and UMAP embeddings, and clustering labels are app specific

However we included the niche labels as: ------

## Pathway activities

_PROGENy_ pathway activity scores can be visualized as a continuous variable if searched with the prefix _progeny_ 

## Deconvolution results

_cell2location_ scores can be visualized as a continuous variable if searched with the prefix _c2l_ 

Deconvolution scores transformed to proportions can be visualized as a continuous variable if searched with the prefix _c2l_props_ 

