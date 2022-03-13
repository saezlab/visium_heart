library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(tidyverse)

sc_data <- loadHDF5SummarizedExperiment("./processed_visium/integration/integrated_slides_sce/")
umap_info <- readRDS("./results/niche_mapping/ct_niches/umap_compositional.rds")
cluster_info <- readRDS("./results/niche_mapping/ct_niches/niche_annotation_ct.rds")
integrated_compositions <- readRDS("./results/niche_mapping/ct_niches/integrated_compositions.rds")
pat_info <- readRDS("./markers/visium_patient_anns_revisions.rds") %>%
  dplyr::select(sample_id, patient_group, patient_id)

# Add UMAP of ILR to the integrated object

meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  mutate(row_id = strsplit(row_id, "_") %>%
           map_chr(., ~.x[[1]])) %>%
  mutate(row_id = paste0(orig.ident, "..", row_id))

colnames(sc_data) <- meta_data$row_id

# Add niche labels

rownames(cluster_info) <- cluster_info$row_id
sc_data$niche <- cluster_info[colnames(sc_data), "ct_niche"]

# Add patient annotations
colData(sc_data) <- cbind(colData(sc_data), integrated_compositions[colnames(sc_data),])

pat_labels <- colData(sc_data) %>% 
  as.data.frame() %>%
  dplyr::select(orig.ident) %>%
  left_join(pat_info, by = c("orig.ident" = "sample_id")) %>%
  dplyr::select(-orig.ident)

colData(sc_data) <- cbind(colData(sc_data), pat_labels)

# Filter slides that are not useful (not niche representation)

niche_props <- colData(sc_data) %>%
  as.data.frame() %>%
  group_by(patient_id, niche)  %>%
  summarise(n_spots = n()) %>%
  dplyr::mutate(n_spots_sample = sum(n_spots)) %>%
  dplyr::mutate(prop_spots = n_spots/n_spots_sample)

niche_distribution <- ggplot(niche_props, aes(x = niche, 
                                              y = patient_id, 
                                              size = prop_spots)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("./results/niche_mapping/ct_niches/niche_distribution_samples.pdf", height = 5, width = 4)

plot(niche_distribution)

dev.off()

high_qc_pats <- niche_props %>%
  group_by(patient_id) %>%
  summarize(max_comp =  max(prop_spots)) %>%
  dplyr::filter(max_comp < 0.8) %>%
  pull(patient_id)

ix <- which(sc_data$patient_id %in% high_qc_pats)

sc_data <- sc_data[, ix]

# Create pseudobulk

pb_data <- scuttle::summarizeAssayByGroup(x = sc_data,
                                          ids = colData(sc_data)[,c("patient_id", "niche")],
                                          statistics = "sum")

saveRDS(pb_data, file = "./processed_visium/integration/pb_nichepat_dat.rds")


