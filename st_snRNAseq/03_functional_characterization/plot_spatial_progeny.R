# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script I plot progeny scores for all slides

path <- "./processed_visium/objects/"
outpath <- "./results/progeny/spatial/"
slide_files <- list.files(path)
slide_names <- gsub("[.]rds", "", slide_files)

param_df <- tibble(slide_name = slide_names,
                   slide_file = paste0(path, slide_files),
                   slide_out = paste0(outpath, slide_files %>% gsub("[.]rds", "_progenyplts.pdf",.)))

plt_spatialprogeny <- function(slide_name, slide_file, slide_out) {
  
  print(slide_name)
  
  visium_slide <- readRDS(slide_file)
  
  DefaultAssay(visium_slide) <- "progeny"
  
  pdf(slide_out, height = 20, width = 16)
  
  prog_plt <- SpatialFeaturePlot(visium_slide, 
                                 features = rownames(visium_slide),
                                 ncol = 4)
  
  plot(prog_plt)
  
  dev.off()
  
}

pwalk(param_df, plt_spatialprogeny)












