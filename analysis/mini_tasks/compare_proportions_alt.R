# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' How correlated are the compositions in visium and rnaSEQ

library(tidyverse)
library(Seurat)

# Get all data files -----------------------------------------------------------------------
slide_files_folder <- "./visium_results_manuscript/processed_visium_revisions/"
slide_files <- list.files(slide_files_folder)
slide_files_full <- paste0(slide_files_folder, 
                           slide_files)
slide_ids <- gsub("[.]rds", "", slide_files)

# Patient annotation -------------------------------------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(Visium, New.Ids) %>%
  dplyr::rename(slide = Visium,
                patient = New.Ids)

# Tibble of all data
# c2l: cell2location
# c2l_states: spotlight

deconv_res <- tibble("slide_path" = slide_files_full,
                     "slide" = slide_ids) %>%
  dplyr::mutate(deconv_mats = map(slide_path, function(x) {
    
    slide_obj <- readRDS(x)
    
    list("cell2location" = as.matrix(slide_obj[["c2l"]]@data) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type),
           "cell2location_states" = as.matrix(slide_obj[["c2l_statesred"]]@data) %>%
             as.data.frame() %>%
             rownames_to_column("cell_type") %>%
             pivot_longer(-cell_type))
    
    
  }))

# Here I keep in long format for further filtering ----------

deconv_res <- deconv_res %>% 
  dplyr::select(slide, deconv_mats) %>%
  dplyr::mutate(cell2location_states = map(deconv_mats, ~ .x[["cell2location_states"]]),
                cell2location = map(deconv_mats, ~ .x[["cell2location"]])) %>%
  dplyr::select(slide, cell2location_states)  %>%
  dplyr::left_join(condition_dictionary)

deconv_res <- deconv_res %>% 
  unnest() %>%
  dplyr::rename(.,c2l_value = value)

# Here we identify spots that contain 0 cells

spot_dict <- deconv_res %>%
  group_by(slide, name) %>%
  summarize(spot_cells =  sum(c2l_value))

deconv_res <- deconv_res %>%
  left_join(spot_dict, by  = c("slide", "name"))

deconv_res <- deconv_res %>%
  dplyr::filter(spot_cells >= 1)

# Here I transform c2l scores to proportions ----------------------------------------
deconv_res <- deconv_res %>%
  group_by(patient, name) %>%
  dplyr::mutate(c2l_prop = c2l_value/(sum(c2l_value))) %>%
  ungroup()

slide_cell_count <- deconv_res %>%
  group_by(patient, cell_type) %>%
  summarise(c2l_cell_prop = sum(c2l_value)) %>%
  ungroup()

all_props <- slide_cell_count %>% 
  group_by(patient) %>%
  mutate(c2l_slide_count = sum(c2l_cell_prop)) %>%
  mutate(c2l_cell_prop = c2l_cell_prop/c2l_slide_count) %>%
  unique()

# Patient annotation -------------------------------------------------------------------------
condition_dictionary <- read.table("./markers/NEW_PatIDs_visium_overview_allsamples.tsv",
                                   sep = "\t", header = T) %>%
  mutate_all(as.character) %>%
  dplyr::select(snRNA, New.Ids) %>%
  dplyr::rename(patient = New.Ids)

# Scell data meta --------------------------------------------------------------------------
scell_meta <- readRDS("visium_results_manuscript/integration/ps_integrated_data_fordeconv.rds")[[1]][["annotations"]] %>%
  dplyr::mutate(deconv_col = gsub("_","-", deconv_col))

# Generate cell type counts -----------------------------------------------------------
state_prop <- scell_meta %>%
  group_by(orig.ident) %>%
  mutate(n_cells = length(orig.ident)) %>%
  ungroup() %>%
  group_by(orig.ident, deconv_col) %>%
  mutate(n_state = length(orig.ident)) %>%
  dplyr::select(orig.ident, 
                n_cells, 
                deconv_col, 
                n_state) %>%
  unique() %>%
  dplyr::rename(snRNA = orig.ident) %>%
  left_join(condition_dictionary)

# Generate cell state proportions --------------------------------------------------------
state_prop <- state_prop %>% 
  ungroup() %>%
  dplyr::mutate(n_state = n_state/n_cells) %>%
  dplyr::select(deconv_col, n_state, patient)

all_props <- all_props %>% 
  left_join(state_prop, by = c("patient", 
                               "cell_type" = "deconv_col"))

all_props_cor <- all_props %>%
  na.omit() %>%
  group_by(patient) %>%
  nest() %>%
  mutate(cstates_cor = map(data, function(x) { 
    
    dat <- na.omit(x)
    broom::tidy(cor.test(log10(dat$c2l_cell_prop), 
                         log10(dat$n_state), 
                         method = "spearman"))
    
    }))

ggplot(all_props, aes(x = log10(c2l_cell_prop), 
                      y = log10(n_state), 
                      color = cell_type)) +
  geom_point() +
  xlab("Visium") +
  ylab("snRNA") +
  facet_wrap(.~patient, nrow =  3)

# Grouping by cell type --------------------------------------------------------

cell_dictionary <- scell_meta %>%
  dplyr::select(cell_type, deconv_col) %>%
  unique() %>%
  dplyr::mutate(cell_type = ifelse(deconv_col == "VSMCs",
                                   "VSMCs", cell_type)) %>%
  dplyr::rename(mjr_ct = cell_type,
                cell_type = deconv_col)

ct_props_grouped <- all_props %>% 
  left_join(cell_dictionary) %>%
  group_by(patient, mjr_ct) %>%
  summarise(mjr_visium_prop = sum(c2l_cell_prop),
            mjr_scell_prop = sum(n_state))

ggplot(ct_props_grouped, aes(x = log10(mjr_visium_prop), 
                        y = log10(mjr_scell_prop), 
                        color = mjr_ct)) +
  geom_point() +
  xlab("Visium") +
  ylab("snRNA") +
  facet_wrap(.~patient, nrow =  3)
  
# Do visium proportions correlate with snRNA per patient?

ct_props_grouped %>%
  ungroup() %>%
  na.omit() %>%
  group_by(patient) %>%
  nest() %>%
  dplyr::mutate(cor_cts = map(data, function(x) {
    
    broom::tidy(cor.test(log10(x$mjr_visium_prop), 
                         log10(x$mjr_scell_prop),
                         method = "spearman"))
    
  })) %>%
  unnest(cor_cts) %>%
  ggplot(aes(x = estimate, y =  -log10(p.value))) +
  geom_point() +
  xlab("Spearman correlation")
    
# Are proportions per cell-type more robust?

ct_props_dat <- all_props %>% 
  left_join(cell_dictionary) %>%
  ungroup() %>%
  group_by(mjr_ct) %>%
  nest() %>%
  mutate(ct_props_plts = map(data, function(x) {
    
    ggplot(x, aes(x = c2l_cell_prop, 
               y = n_state, 
               color = cell_type)) +
      geom_point() +
      xlab("Visium") +
      ylab("snRNA") +
      facet_wrap(.~patient, nrow =  3, 
                 scales = "free")
  
  })) %>%
  mutate(ct_props_cor = map(data, function(x) { 
    
    dat <- na.omit(x)
    broom::tidy(cor.test(dat$c2l_cell_prop, dat$n_state, method = "spearman"))
    
    
    }))

pdf("visium_results_manuscript/ct_data//evaluate_deconv_props_states_red.pdf")

walk(ct_props_dat$ct_props_plts, print)

dev.off()

ct_props_dat %>%
  dplyr::filter(mjr_ct == "macrophages") %>%
  unnest(data) %>%
  ggplot(aes(x = c2l_cell_prop, y = n_state)) +
  geom_point()
  









