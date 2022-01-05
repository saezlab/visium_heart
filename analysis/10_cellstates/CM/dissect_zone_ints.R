# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Get only BZ

slide_annotation <- read_csv("./markers/visium_patient_anns_revisions.csv")
misty_out_folder <- "./results/state_structure/CM_ct/"

plot_hmaps <- function(misty_res_slide) {
  
  mistyR::plot_improvement_stats(misty_res_slide)
  mistyR::plot_view_contributions(misty_res_slide)
  
  view_names <- misty_res_slide$importances$view %>% unique()
  walk(view_names, function(v) {
    mistyR::plot_interaction_heatmap(misty_res_slide, v, cutoff = 0)
  })
  
}

BZ <- slide_annotation %>%
  dplyr::filter(major_labl == "BZ") %>%
  pull(sample_id) %>%
  paste0(misty_out_folder, "mstate_", .) %>%
  mistyR::collect_results(.)

pdf("./results/state_structure/CM_ct/BZ_summary.pdf", height = 5, width = 5)

plot_hmaps(BZ)

dev.off()


RZ <- slide_annotation %>%
  dplyr::filter(major_labl == "RZ") %>%
  pull(sample_id) %>%
  paste0(misty_out_folder, "mstate_", .) %>%
  mistyR::collect_results(.)

pdf("./results/state_structure/CM_ct/RZ_summary.pdf", height = 5, width = 5)

plot_hmaps(RZ)

dev.off()

CTRL <- slide_annotation %>%
  dplyr::filter(major_labl == "CTRL") %>%
  pull(sample_id) %>%
  paste0(misty_out_folder, "mstate_", .) %>%
  mistyR::collect_results(.)

pdf("./results/state_structure/CM_ct/CTRL_summary.pdf", height = 5, width = 5)

plot_hmaps(CTRL)

dev.off()







