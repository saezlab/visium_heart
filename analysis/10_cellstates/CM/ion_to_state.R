# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we check to which states
#' the ordered probes belong

library(tidyverse)
library(cowplot)
library(Seurat)
library(ggpubr)
source("./analysis/utils/funcomics.R")
source("./analysis/utils/misty_pipeline.R")

sample_dict <- read_csv("./markers/visium_patient_anns_revisions.csv")

# Generate the dictionary of channel sets ---------------------------------------------------

probes <- read_csv("./markers/ion_channels.csv",col_names = T) %>%
  dplyr::rename("channel" = `transmural ion channel genes`) %>%
  unique() %>%
  dplyr::pull(channel)

markers <- list("REACTOME_ION_CHANNEL_TRANSPORT" = readRDS("./markers/Genesets_Dec19.rds")$MSIGDB_CANONICAL[["REACTOME_ION_CHANNEL_TRANSPORT"]],
                "TRANSMURAL_ION_CHANNELS" = probes)


# Get individual slide info ---------------------------------------------
visium_folder = "./processed_visium/objects/"
out_folder = "./results/ion_plots/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)
visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    sample = visium_samples,
                    out_file = paste0(out_folder, 
                                      sample,"_ionchannel.pdf"))

gsets <- enframe(markers,name = "source", value = "target") %>% unnest() %>% mutate(mor = 1, likelihood = 1)


state_ions <- map(set_names(visium_df$visium_file, visium_samples), 
    function(visium_file) { 
  
  print(visium_file)
  
  visium_slide <- readRDS(visium_file) %>%
    positive_states(., assay = "cell_states") %>%
    filter_states(slide = .,
                  by_prop = F,
                  prop_thrsh = 0.10)
  
  visium_slide <- visium_slide %>%
    get_wmean_score(visium_slide = .,
                    network = gsets,
                    assay = "SCT",
                    module_name = "ion_channels")
  
  DefaultAssay(visium_slide) <- "Spatial"
  
  rbind(GetAssayData(visium_slide, assay = "ion_channels"), 
        GetAssayData(visium_slide, assay = "cell_states_pos"))
  
})


state_ions_cors <- map(state_ions, function(dat) {
  
  filt_data <- dat %>%
    t() %>%
    as.data.frame() %>%
    select_at(c("REACTOME-ION-CHANNEL-TRANSPORT", "TRANSMURAL-ION-CHANNELS", 
                "CM-damaged-CM","CM-healthy-CM", "CM-intermediate-CM")) %>%
    dplyr::filter(`CM-damaged-CM` > 0 |
                    `CM-healthy-CM` > 0 |
                    `CM-intermediate-CM` > 0)
  
  ions <- colnames(filt_data)[grepl(c("ION"), colnames(filt_data))]
  states <- colnames(filt_data)[grepl(c("CM"), colnames(filt_data))]
  
  cor_res <- crossing(ions_names = ions, 
                      states_names = states) %>%
    mutate(sprmn_res = map2(ions_names, states_names, function(i, s) {
      
      cor.test(filt_data[,i], filt_data[,s], method = "spearman") %>%
        broom::tidy()
      
    })) %>%
    unnest()
  
})

# Compare the correlations between ion channels mean scores
# and states for each slide

state_order <- c("CM-healthy-CM", "CM-intermediate-CM", "CM-damaged-CM")

my_comparisons <- list( c("CM-healthy-CM", "CM-intermediate-CM"), 
                        c("CM-healthy-CM", "CM-damaged-CM"), 
                        c("CM-intermediate-CM", "CM-damaged-CM") )


corr_plts <- enframe(state_ions_cors) %>%
  unnest() %>%
  dplyr::mutate(ions_names = ifelse(ions_names == "REACTOME-ION-CHANNEL-TRANSPORT",
                                    "ion-channel-transport",
                                    "transmural-ion_channels"),
                states_names = factor(states_names, levels = state_order)) %>%
  
  ggplot(aes(x = states_names, y = estimate)) +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
        ) +
  stat_compare_means(comparisons = my_comparisons, size = 2.1) + # Add pairwise comparisons p-value
  facet_wrap(. ~ ions_names) +
  ylab("Spearman correlation") +
  xlab("")

pdf("./results/ion_plots/cor_ions_states.pdf", height = 5, width = 4)

plot(corr_plts)

dev.off()

write.csv(enframe(state_ions_cors) %>%
            unnest(), "./results/ion_plots/cor_ions_states.csv")

saveRDS(state_ions, file = "./results/ion_plots/raw_ions_states.rds")








