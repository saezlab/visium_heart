# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will add cell-death signatures to the slides
#' 
#' Our assumption is that in all the IZ slides we will have a greater expression

library(Seurat)
library(tidyverse)
source("./analysis/utils/funcomics.R")

# First, prove that slides with more loss of spots
# have less cells

snrna_qc <- read_csv("./processed_snrnaseq/initial_qc/all_qcs_after_filtering_snrna.csv")
spatial_qc <- read_csv("./processed_visium/initial_qc/all_qcs_after_filtering.csv")

pat_meta <- left_join(snrna_qc, spatial_qc) %>%
  dplyr::mutate(perc_spot_loss = 1 - (nspots/n_spots_under_tissue))

patient_ann  <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  left_join(pat_meta)

pat_meta <- pat_meta %>% 
  dplyr::mutate(loss = ifelse(perc_spot_loss >= 0.025,
                              TRUE, FALSE))

pdf("results/cell_death/spot_loss_ncells.pdf", height = 4, width = 4)
print(ggplot(pat_meta, aes(x = log10(ncells), y = perc_spot_loss)) + geom_point())

print(ggplot(pat_meta, aes(y = log10(ncells), x = loss, color = loss)) + geom_boxplot() +
  geom_jitter())
dev.off()

# First find gene sets of interest
gsets <- readRDS("./markers/Genesets_Dec19.rds")
gsets <- c(gsets$MSIGDB_HMARKS, gsets$MSIGDB_CANONICAL, gsets$MSIGDB_GO_BIOLPROC)
  
apoptosis <- gsets[grepl("apoptosis", 
                         names(gsets),
                         ignore.case = T)]
necrosis <- gsets[grepl("necrosis", 
                        names(gsets),
                        ignore.case = T)]

death <- gsets[grepl("death", 
                     names(gsets),
                     ignore.case = T)]

death_gsets <- c(apoptosis, necrosis, death)

# Add these genesets into filtered slides...

death_gsets <- enframe(death_gsets) %>%
  unnest() %>%
  mutate(mor = 1,
         likelihood = 1) %>%
  dplyr::rename("source" = name,
                "target" = value) %>%
  unique()


estimate_death_scores <- function(slide_file) {
  print(slide_file)
  
  visium_slide <- readRDS(slide_file)
  
  visium_slide <- get_wmean_score(visium_slide = visium_slide,
                                  network = death_gsets,
                                  assay = "Spatial",
                                  module_name = "death_gsets")
  
  GetAssayData(visium_slide, assay = "death_gsets")
  
}

visium_folder <- "/Volumes/RicoData/mi_atlas_extradata/processed_visium_unfiltered/objects/"
sample_file <- list.files(visium_folder)

param_df <- tibble(sample_id = gsub("[.]rds", "", 
                                    sample_file),
                   slide_file = paste0(visium_folder,
                                       sample_file))

cell_death_scores <- param_df %>%
  mutate(death_scores = map(slide_file, estimate_death_scores))

saveRDS(cell_death_scores, file = "./results/cell_death/cell_death_scores.rds")

# We can try the simplest which is to take the mean 
# If more spots contain these death regions then they will drag the mean

n_cells_plt <- patient_ann %>%
  dplyr::select(patient_id, patient_group, ncells) %>%
  unique() %>%
  ggplot(aes(x = patient_group, y = ncells, color = patient_group)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitter()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  ggtitle("number of cells")

pdf("./results/cell_death/cell_numbers.pdf", height = 4, width = 3)

plot(n_cells_plt)

dev.off()

prop_scores <- map(set_names(cell_death_scores$death_scores, 
              cell_death_scores$sample_id), function(x) {
                pos_set <- (x > 1)
                m_s <- rowMeans(pos_set, na.rm = T)
                m_s <- m_s/ncol(x)
       tibble(gset = names(m_s), score = m_s)
     }) %>%
  enframe() %>%
  unnest() %>%
  left_join(patient_ann %>% dplyr::select(sample_id, patient_group) %>% unique(),
            by = c("name" = "sample_id")) %>%
  group_by(gset) %>%
  nest()

mean_scores <- map(set_names(cell_death_scores$death_scores, 
                             cell_death_scores$sample_id), function(x) {
                               x[is.infinite(x)] <- NA
                               m_s <- rowMeans(x, na.rm = T)
                               tibble(gset = names(m_s), score = m_s)
                             }) %>%
  enframe() %>%
  unnest() %>%
  left_join(patient_ann %>% dplyr::select(sample_id, patient_group) %>% unique(),
            by = c("name" = "sample_id")) %>%
  group_by(gset) %>%
  nest()

prop_inf <- map(set_names(cell_death_scores$death_scores, 
                          cell_death_scores$sample_id), function(x) {
                            pos_set <- is.infinite(x)
                            m_s <- rowSums(pos_set, na.rm = T)
                            m_s <- m_s/ncol(x)
                            tibble(gset = names(m_s), score = m_s)
                          }) %>%
  enframe() %>%
  unnest() %>%
  left_join(patient_ann %>% dplyr::select(sample_id, patient_group) %>% unique(),
            by = c("name" = "sample_id")) %>%
  group_by(gset) %>%
  nest()


pdf("./results/cell_death/functional_meanscores.pdf", height = 4, width = 3)

walk2(mean_scores$gset, mean_scores$data, function(g,d) {
  
  b_plt <- d %>%
    ggplot(aes(x = patient_group, y = score, color = patient_group)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, position = position_jitter()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    ggtitle(g)
  
  plot(b_plt)
  
})

dev.off()

pdf("./results/cell_death/functional_propscores.pdf", height = 4, width = 3)

walk2(prop_scores$gset, prop_scores$data, function(g,d) {
  
  b_plt <- d %>%
    ggplot(aes(x = patient_group, y = score, color = patient_group)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(pch = 21, position = position_jitter()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") + 
    ggtitle(g)
  
  plot(b_plt)
  
})

dev.off()

pdf("./results/cell_death/infinite_props.pdf", height = 4, width = 3)

walk2(prop_inf$gset, prop_inf$data, function(g,d) {
  
  b_plt <- d %>%
    ggplot(aes(x = patient_group, y = score, color = patient_group)) +
    geom_boxplot() +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") + 
    ggtitle(g)
  
  plot(b_plt)
  
})

dev.off()
