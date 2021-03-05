# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Copying cell2location models to simpler paths
library(tidyverse)

root_path <- "./results/deconvolution/c2l/"
c2l_files <- c("LocationModelLinearDependentW_1experiments_9clusters_3011locations_3426genes_157782/",
               "LocationModelLinearDependentW_1experiments_9clusters_3086locations_3426genes_157785/",
               "LocationModelLinearDependentW_1experiments_9clusters_3377locations_3426genes_157777/",
               "LocationModelLinearDependentW_1experiments_9clusters_3774locations_3426genes_157779/",
               "LocationModelLinearDependentW_1experiments_9clusters_4269locations_3426genes_157771/",
               "LocationModelLinearDependentW_1experiments_9clusters_4270locations_3426genes_157775/",
               "LocationModelLinearDependentW_1experiments_9clusters_4548locations_3426genes_157772/",
               "LocationModelLinearDependentW_1experiments_9clusters_4661locations_3426genes_157781/")
samples <- gsub("/","",map_chr(strsplit(c2l_files,"_"), last))
plots_file <- "plots"
density <- "W_cell_density.csv"  
density_q <- "W_cell_density_q05.csv"


# Copy the plots folder of each model to location_models/plots
walk2(c2l_files, samples, function(c2l_f, s) {
  
  print(c2l_f)
  
  # First copy as it is
  system(paste0("cp -R ", 
                root_path, 
                c2l_f, 
                plots_file,
                " ./results/deconvolution/c2l/location_models/plots/"
                ))
  
  # Then rename
  
  paste0(" ./results/deconvolution/c2l/location_models/plots/",
         plots_file,
         " ./results/deconvolution/c2l/location_models/plots/",
         plots_file,
         "_",
         s)
  
  system(paste0("mv",
                " ./results/deconvolution/c2l/location_models/plots/",
                plots_file,
                " ./results/deconvolution/c2l/location_models/plots/",
                plots_file,
                "_",
                s))
  
})

# Copy the density file of each model to location_models/density_tables

walk2(c2l_files, samples, function(c2l_f, s) {
  
  print(c2l_f)
  
  # First copy
  
  system(paste0("cp ", 
                root_path, 
                c2l_f, 
                density,
                " ./results/deconvolution/c2l/location_models/density_tables/"))
  
  # Then move
  
  system(paste0("mv",
                " ./results/deconvolution/c2l/location_models/density_tables/",
                density,
                " ./results/deconvolution/c2l/location_models/density_tables/",
                s,
                "_",
                density))
  
})

# Copy the density file of each model to location_models/density_tables

walk2(c2l_files, samples, function(c2l_f, s) {
  
  print(c2l_f)
  
  system(paste0("cp ", 
                root_path, 
                c2l_f, 
                density_q,
                " ./results/deconvolution/c2l/location_models/density_tables/"))
  
  system(paste0("mv",
                " ./results/deconvolution/c2l/location_models/density_tables/",
                density_q,
                " ./results/deconvolution/c2l/location_models/density_tables/",
                s,
                "_",
                density_q))
  
})
