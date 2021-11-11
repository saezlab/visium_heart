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

sc_data <- loadHDF5SummarizedExperiment("./results/ct_data/CM/CM_states_sce/")

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

pdf("./results/ct_data/CM/qc_states.pdf")
plot(qc_panel)
dev.off()

# Parameters ---------------------------------------------------------------------

# Global ----------------------------------
patient_annotations <- readRDS("./markers/snrna_patient_anns_revisions.rds")
opt_state <- "opt_state" # This is always the same
group_alias <- "state"

# Load SCE object
ct_folder <- "./results/ct_data/"
ct_alias <- "CM"
exception <- c(1, 4,7) # Here you define the cluster to exclude 
perc_thrsh <- 0.1 #Minimum percentage of expression of a gene within a state to be considered
n_samples_filt <- 5 # Minimum number of samples in each state

all_panels <- summarize_state(ct_folder, ct_alias, 
                              exception, perc_thrsh, 
                              n_samples_filt, background_exclude = "CM")

# Annotate this with HCA

mi_markers <- readRDS("./results/ct_data/CM/state_genelist.rds")

hca_markers <- read_csv("./ext_data/hca_cardiac_mrkrs.csv") %>%
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

GSE_analysis = function(geneList,Annotation_DB){
  
  geneList = geneList[geneList %in% unique(unlist(Annotation_DB))]
  
  ResultsDF = matrix(0,nrow = length(Annotation_DB),ncol = 5)
  rownames(ResultsDF) = names(Annotation_DB)
  colnames(ResultsDF) = c("GenesInPathway","GenesInList","GeneNames","p_value","corr_p_value")
  
  DB_genecontent = length(unique(unlist(Annotation_DB)))
  
  GenesDB = DB_genecontent 
  SelectedGenes = length(geneList)
  
  for(gset in rownames(ResultsDF)){
    GP = length(Annotation_DB[[gset]])
    GL = length(intersect(Annotation_DB[[gset]],geneList))
    
    ResultsDF[gset,"GenesInList"] = GL
    ResultsDF[gset,"GenesInPathway"] = GP
    ResultsDF[gset,"GeneNames"] = paste(intersect(Annotation_DB[[gset]],geneList),collapse = ",")
    ResultsDF[gset,"p_value"] = phyper(q=GL - 1, m=GP, n=GenesDB-GP, k=SelectedGenes, lower.tail = FALSE, log.p = FALSE)
  }
  
  ResultsDF[,"corr_p_value"] = p.adjust(ResultsDF[,"p_value"],method = "BH")
  ResultsDF = data.frame(ResultsDF,stringsAsFactors = F)
  ResultsDF = ResultsDF[order(ResultsDF[,"p_value"]),]
  
  ResultsDF = ResultsDF %>% 
    rownames_to_column("gset") %>% 
    mutate_at(c("GenesInPathway","GenesInList",
                "p_value","corr_p_value"), 
              as.numeric) %>% 
    dplyr::arrange(corr_p_value,GenesInList)
  
  return(ResultsDF)
  
}

label_mapping <- map(mi_markers, GSE_analysis, Annotation_DB = hca_markers) %>% 
  enframe() %>%
  unnest() %>%
  dplyr::select(name, gset, corr_p_value)

label_mapping_plt <- ggplot(label_mapping, aes(x = name, 
                                               y = gset, 
                                               size = -log10(corr_p_value))) +
  geom_point()


pdf("./results/ct_data/Endo/label_mapping_plt.pdf")

plot(label_mapping_plt)

dev.off()


# Re do it again

source("./analysis/10_cellstates/get_state_summary.R")

# Data evaluation ----------------------------------------------------------------
# Here we do human cutoffs to define which cells to keep for state definition ----

sc_data <- loadHDF5SummarizedExperiment("./results/ct_data_filt/CM/CM_states_sce/")

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

pdf("./results/ct_data/CM/qc_states.pdf")
plot(qc_panel)
dev.off()

# Parameters ---------------------------------------------------------------------

# Global ----------------------------------
patient_annotations <- readRDS("./markers/snrna_patient_anns_revisions.rds")
opt_state <- "opt_state" # This is always the same
group_alias <- "state"

# Load SCE object
ct_folder <- "./results/ct_data_filt/"
ct_alias <- "CM"
exception <- c(3,4,5,6) # Here you define the cluster to exclude 
perc_thrsh <- 0.1 #Minimum percentage of expression of a gene within a state to be considered
n_samples_filt <- 5 # Minimum number of samples in each state

all_panels <- summarize_state(ct_folder, ct_alias, 
                              exception, perc_thrsh, 
                              n_samples_filt, background_exclude = "CM")

# Annotate this with HCA

mi_markers <- readRDS("./results/ct_data_filt/CM/state_genelist.rds")

hca_markers <- read_csv("./ext_data/hca_cardiac_mrkrs.csv") %>%
  dplyr::filter(logfoldchanges > 0, pvals_adj < 0.001) %>%
  dplyr::arrange(group, -logfoldchanges) %>%
  dplyr::group_by(group) %>%
  dplyr::slice(1:200) %>%
  dplyr::select(group, names) %>%
  nest() %>%
  dplyr::mutate(data = map(data, ~ .x[[1]])) %>%
  deframe()
