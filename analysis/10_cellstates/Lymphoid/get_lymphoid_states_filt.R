# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform a whole cell state characterization
#' 
#' 1) We identify genes that are representative for each state (40% cutoff)
#' 2) We check proportions of states in patients and keep states that are representative of a population (larger than 1% in 5 patients)
#' 3) After state filtering we perform linear mixed models to identify genes of interest
#' 4) We used those t-values for pathway activity estimation and tf activity estimation

source("./analysis/10_cellstates/get_state_summary.R")

# Data evaluation ----------------------------------------------------------------
# Here we do human cutoffs to define which cells to keep for state definition ----

sc_data <- loadHDF5SummarizedExperiment("./results/ct_data_filt/Lymphoid/Lymphoid_states_sce/")

meta_data <- colData(sc_data) %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  dplyr::mutate(cell_id = strsplit(cell_id, "_") %>%
                  map_chr(., ~ .x[[1]]))

doublet_score <- meta_data %>%
  ggplot(aes(y = log(doublet_score), x = opt_state)) +
  geom_boxplot()

features <- meta_data %>%
  ggplot(aes(y = nFeature_RNA, x = opt_state)) +
  geom_boxplot()

counts <- meta_data %>%
  ggplot(aes(y = nCount_RNA, x = opt_state)) +
  geom_boxplot()

mt <- meta_data %>%
  ggplot(aes(y = log(percent.mt), x = opt_state)) +
  geom_boxplot()

qc_panel <- cowplot::plot_grid(doublet_score, features,
                               counts, mt, ncol = 2, align = "hv")

pdf("./results/ct_data_filt/Lymphoid/qc_states.pdf")
plot(qc_panel)
dev.off()

#Parameters ---------------------------------------------------------------------

# Global ----------------------------------
patient_annotations <- readRDS("./markers/snrna_patient_anns_revisions.rds")
opt_state <- "opt_state" # This is always the same
group_alias <- "state"

# Load SCE object
ct_folder <- "./results/ct_data_filt/"
ct_alias <- "Lymphoid"
exception <- NA # Here you define the cluster to exclude 
perc_thrsh <- 0.1 #Minimum percentage of expression of a gene within a state to be considered
n_samples_filt <- 5 # Minimum number of samples in each state

all_panels <- summarize_state(ct_folder, ct_alias, 
                              exception, perc_thrsh, 
                              n_samples_filt, background_exclude = "Lymphoid")

# Run 

mi_markers <- readRDS("./results/ct_data_filt/Lymphoid/state_genelist.rds")

hca_markers <- read_csv("./ext_data/hca_immune_mrkrs.csv") %>%
  dplyr::filter(logfoldchanges > 0, pvals_adj < 0.001) %>%
  dplyr::arrange(group, -logfoldchanges) %>%
  dplyr::group_by(group) %>%
  dplyr::slice(1:200) %>%
  dplyr::select(group, names) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()

# Now enrich with hypergeometric tests

# Function to do enrichment -----------------------------------------------------------------

label_mapping <- map(mi_markers, GSE_analysis, Annotation_DB = hca_markers) %>% 
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, corr_p_value) %>%
  dplyr::filter(corr_p_value < 0.05)


label_mapping_plt <- ggplot(label_mapping, aes(x = name, 
                                               y = gset, 
                                               size = -log10(corr_p_value))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


pdf("./results/ct_data_filt/Lymphoid/label_mapping_plt.pdf", height = 4, width = 5)

plot(label_mapping_plt)

dev.off()
