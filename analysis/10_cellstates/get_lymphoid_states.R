# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we perform a whole cell state characterization
#' 
#' 1) We identify genes that are representative for each state (40% cutoff)
#' 2) We check proportions of states in patients and keep states that are representative of a population (larger than 1% in 5 patients)
#' 3) After state filtering we perform linear mixed models to identify genes of interest
#' 4) We used those t-values for pathway activity estimation and tf activity estimation

source("./analysis/10_cellstates/get_state_summary.R")

# Parameters ---------------------------------------------------------------------

# Global ----------------------------------
patient_annotations <- readRDS("./markers/snrna_patient_anns_revisions.rds")
opt_state <- "opt_state" # This is always the same
group_alias <- "state"

# Load SCE object
ct_folder <- "./results/ct_data/"
ct_alias <- "Lymphoid"
exception <- c(4, 5) # Here you define the cluster to exclude 
perc_thrsh <- 0.2 #Minimum percentage of expression within a state to be considered
n_samples_filt <- 5 # Minimum number of samples in each state

all_panels <- summarize_state(ct_folder, ct_alias, exception, perc_thrsh, n_samples_filt)



