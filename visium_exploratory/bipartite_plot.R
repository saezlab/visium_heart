library(tidyverse)
library(bipartite)
library(RColorBrewer)

bipartite_plot <- function(result, threshold){

  mat <- pivot_wider(result[["para_importance"]], names_from = Predicted,  
                     values_from=Importance, values_fill = 0) %>%
    column_to_rownames("Predictor")
  
  # if shared the marker will be assigned the first available ct
  cts <- left_join(unique(result[["para_importance"]] %>% select(Predicted)),
                   result[["marker_dictionary"]], by = c("Predicted" = "gene"))
  
  mat <- mat[,unique(cts %>% pull(Predicted))]
  mat[mat <= threshold] <- 0
  
  colorids <- as.numeric(as.factor(cts$ct))
  colors <- brewer.pal(length(unique(cts$ct)),"Set1")
  
  plotweb(mat, col.high = colors[colorids])
  legend("right", levels(as.factor(cts$ct)), fill = colors)
}
