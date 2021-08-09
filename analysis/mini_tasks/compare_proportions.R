# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' How different are cell_proportions
#' 
#' 

library(tidyverse)

# Patient annotation -------------------------------------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character)

# Scell data --------------------------------------------------------------------------
scell_meta <- readRDS("visium_results_manuscript/integration/integrated_wstates_meta.rds")

state_prop <- scell_meta %>%
  group_by(orig.ident) %>%
  mutate(n_cells = length(orig.ident)) %>%
  ungroup() %>%
  group_by(orig.ident, deconv_col) %>%
  mutate(n_state = length(deconv_col)) %>%
  dplyr::select(orig.ident, n_cells, deconv_col, n_state) %>%
  unique() %>%
  rename(snRNA = orig.ident) %>%
  left_join(condition_dictionary)

prop_overview <- ggplot(state_prop, 
                        aes(fill = deconv_col, 
                            x = n_state, 
                            y = New.Ids)) + 
  geom_bar(position="fill", 
           stat="identity") +
  theme(legend.position = "bottom") +
  ylab("") + xlab("proportions")   


# Read deconv matrices --------------------
deconv_mats_folder <- "./visium_results_manuscript/deconvolution/spotlight/"
deconv_mats <- list.files(deconv_mats_folder)
deconv_mats <- deconv_mats[grepl(".rds",deconv_mats)]
deconv_mats_ids <- gsub("_deconvmat[.]rds", "", deconv_mats)
deconv_mats <- paste0(deconv_mats_folder, deconv_mats)
visium_prop <- map(set_names(deconv_mats, deconv_mats_ids), readRDS)

visium_prop <- map(visium_prop, function(x) { 
  total <- nrow(x)
  props <- colSums(x)/total
  propr <- tibble(deconv_col = names(props),
                  n_state = props)
  return(propr)
})

visium_prop <- enframe(visium_prop, "Visium") %>%
  unnest() %>%
  left_join(condition_dictionary, by = "Visium") %>%
  mutate(deconv_col = gsub("state-","",deconv_col))


prop_overview <- ggplot(visium_prop, 
                        aes(fill = deconv_col, 
                            x = n_state, 
                            y = New.Ids)) + 
  geom_bar(position="fill", 
           stat="identity") +
  theme(legend.position = "bottom") +
  ylab("") + xlab("proportions")   

# Are proportions similar?

state_prop <- state_prop %>% 
  ungroup() %>%
  dplyr::mutate(n_state = n_state/n_cells) %>%
  dplyr::select(deconv_col, n_state, New.Ids)

visium_prop <- dplyr::select(visium_prop, deconv_col, n_state, New.Ids)

test <- left_join(visium_prop, state_prop, by =c("New.Ids", "deconv_col"))


ggplot(test, aes(x = n_state.x, y = n_state.y, color = deconv_col)) +
  geom_point() +
  xlab("Visium") +
  ylab("snRNA") +
  facet_wrap(.~New.Ids, nrow =  3)
