# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we estimate cell compositions from ATAC
library(tidyverse)

# ATAC data meta --------------------------------------------------------------------------
atac_ann <- readRDS("./markers/atac_patient_anns_revisions.rds")

atac_meta <- read.csv("./processed_atac/metadata.csv") %>%
  dplyr::select(X, orig.ident, cell_type)  %>%
  left_join(atac_ann, by = c("orig.ident" = "sample_id"))# Then fix the annotations ussing the updated table


# Generate cell type counts -----------------------------------------------------------
atac_props <- atac_meta %>%
  group_by(patient_id, cell_type) %>%
  summarize(atac_n_cells = length(patient_id)) %>%
  mutate(atac_prop_cells = atac_n_cells/sum(atac_n_cells)) %>%
  dplyr::mutate(cell_type = gsub("-","_", cell_type)) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "Pericyte", "PC", cell_type)) %>%
  dplyr::mutate(cell_type = ifelse(cell_type == "neuronal", "Neuronal", cell_type)) 

write.table(atac_props, file = "./results/compositions/atac_compositions.txt", col.names = T, row.names = F, quote = F, sep = "\t")

atac_props$keep <- ifelse((atac_props$atac_n_cells > 10), 
                         TRUE, FALSE)

cols <- c("TRUE" = "grey",
          "FALSE" = "black")

# First check the number of cells per patient
patient_summary_plt_flx <- ggplot(atac_props, 
                                  aes(y = patient_id, x = cell_type, fill = keep)) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background= element_rect(fill="black", colour="black")) +
  scale_fill_manual(values = cols,
                    na.value = 'black') +
  coord_equal()

pdf("./results/cell_markers/cell_type_recovery_atac.pdf", height = 5, width = 4)

plot(patient_summary_plt_flx)

dev.off()
