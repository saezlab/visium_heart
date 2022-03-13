# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Utilities to generate ROI based on quadrants


get_quadrant <- function(visium_slide, col_divisions = 4, row_divisions = 4, coord = c(1,1)) {
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium_slide,
                                   cols = c("row", "col"), scale = NULL)
  
  row_min <- min(geometry$row)
  row_max <- max(geometry$row)
  
  row_length <- row_max - row_min
  
  col_min <- min(geometry$col)
  col_max <- max(geometry$col)
  
  col_length <- col_max - col_min
  
  # Make cut-offs
  
  row_cut <- ceiling(row_length/row_divisions)
  col_cut <- ceiling(col_length/col_divisions)
  
  # Now get conditionals
  
  row_max_cut <- (row_min + (row_cut * coord[1]))
  row_min_cut <- (row_min + (row_cut * (coord[1]-1)))
  
  col_max_cut <- (col_min + (col_cut * coord[2]))
  col_min_cut <- (col_min + (col_cut * (coord[2]-1)))
  
  # Filter cells by geometry
  
  ids <- geometry %>%
    as.data.frame() %>%
    rownames_to_column("spot_id") %>%
    filter((row >= row_min_cut) & (row <= row_max_cut),
           (col >= col_min_cut) & (col <= col_max_cut)) %>%
    pull(spot_id)
  
  #return cell_ids
  
  return(ids)
}




