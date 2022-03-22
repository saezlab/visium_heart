# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Utilities for dimensionality reduction


#' Get importances of specific PCs
#' @param prcomp_obj: an object from a prcomp run
#' @param pc: character, it is usually in the form "PC#"
#' @param ntop: numeric. Number of top features to consider based in abs(value)

get_pc_importances <- function(prcomp_obj, pc = "PC1", ntop = 20) {
  
  prcomp_obj$rotation %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    pivot_longer(-term) %>%
    arrange(name, -abs(value)) %>%
    dplyr::filter(name == pc) %>%
    slice(1:ntop) %>%
    arrange(value)
}

#' Create a ggbarplot of importances
#' @param pc_importances: an output object from get_pc_importances
plot_pc_importances <- function(pc_importances) {
  
  pc_importances_order <- pc_importances %>%
    pull(term)
  
  ggplot(pc_importances,
         aes(x = factor(term, levels = pc_importances_order), 
             y = value)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5))  +
    xlab("") + ylab("Importance")
}















