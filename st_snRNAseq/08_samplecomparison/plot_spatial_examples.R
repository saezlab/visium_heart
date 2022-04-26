# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we show the examples of specific interactions

library(tidyverse)
library(ggpubr)
library(Seurat)
library(viridis)

annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

sample_importances_filt <- read.csv("./results/sample_comparison/spatial/sample_importances_filt.csv") %>%
  left_join(annotation_names) %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))) 



plot_importance_boxes <- function(sample_importances_filt, view_sel = "intra", cell_pair) {
  
  my_comparisons <- list( c("myogenic-enriched", "ischemic-enriched"), 
                          c("myogenic-enriched", "fibrotic-enriched"), 
                          c("ischemic-enriched", "fibrotic-enriched") )
  
  
  red_data <- sample_importances_filt %>%
    dplyr::select(name, Predictor, 
                  Target, Importance, 
                  patient_group, view) %>%
    dplyr::filter((Predictor == cell_pair[1] &
                     Target == cell_pair[2]) | 
                    (Predictor == cell_pair[2] &
                       Target == cell_pair[1]),
                  view == view_sel) %>%
  group_by(Predictor, Target) %>%
    nest() %>%
    dplyr::mutate(bplot = map(data, function(dat) {
      
      pw_data <- compare_means(Importance ~ patient_group, 
                               comparisons = my_comparisons, 
                               data = dat) %>%
        dplyr::select(group1, group2, p , p.adj)
      
      
      max_val <- max(dat$Importance) + 0.05
      
      ggplot(dat,
             aes(x = patient_group, color = patient_group, y = Importance)) +
        geom_boxplot() +
        geom_point() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              legend.position = "none",
              axis.text = element_text(size = 11),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        ggpubr::stat_pvalue_manual(pw_data, label = "p.adj", 
                                   y.position = max_val, 
                                   step.increase = 0.1,
                                   tip.length = 0.01,size = 3) +
        #stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
        ggtitle(view_sel) 
    }))
  
    
  return(red_data)
  
}

show_useful_slides <- function(imp_res, pname, tname, pg_name) {
  
  imp_res %>%
    dplyr::select(-bplot) %>%
    unnest() %>%
    filter(patient_group == pg_name) %>%
    arrange(-Importance) %>%
    dplyr::filter(Predictor == pname, Target == tname) 
  
}

plot_visium_panels <- function(plot_data, assay = "c2l", out_folder) {
  
  dir.create(out_folder)
  
  param_df <- plot_data %>%
    mutate(visium_path = paste0("./processed_visium/misty_red_objects/", 
                                name, 
                                "_mistyassays.rds"),
           out_pdf = paste0(out_folder, 
                            "/",
                            name, 
                            "_",
                            patient_group,
                            ".pdf"))
    
  plt_util <- function(Predictor, Target, visium_path, out_pdf) {
    
    print(visium_path)
    
    visium_slide <- readRDS(visium_path)
    
    DefaultAssay(visium_slide) <- assay
    
    plts <- SpatialFeaturePlot(visium_slide, features = c(Predictor, Target),
                               combine = F,max.cutoff = "q99",stroke = 0)
    
    plts <- map(plts, function(spec_f) {
      
      spec_f +
        scale_fill_viridis(option = "A")
  
    })
    
    plts_panel <- cowplot::plot_grid(plotlist = plts, ncol = 2)
    
    pdf(file = out_pdf, height = 4, width = 8)
    
    plot(plts_panel)
    
    dev.off()
    
    return(NULL)
    
  }
  
  pwalk(param_df %>%
          dplyr::select(Predictor, Target, visium_path, out_pdf),
        plt_util)
  
}


# PC to Fibroblast story

# PC to  CM

cm_pc <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                               view_sel = "intra",
                               cell_pair = c("CM", "PC"))

pdf("./results/sample_comparison/spatial/PC_interactions/PCpred_CMtar_intra_box.pdf", height = 4.5, width = 2)

plot(cm_pc$bplot[[2]])

write_csv(cm_pc$data[[2]], 
          file = "./results/sample_comparison/spatial/PC_interactions/PCpred_CMtar_intra_box.csv")

dev.off()


cm_pc_plt_dat <- bind_rows(show_useful_slides(imp_res = cm_pc,
                                              pname = "PC",
                                              tname = "CM",
                                              pg_name = "fibrotic-enriched"),
                           show_useful_slides(imp_res = cm_pc,
                                              pname = "PC",
                                              tname = "CM",
                                              pg_name = "myogenic-enriched") %>%
                             tail())


plot_visium_panels(plot_data = cm_pc_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/PC_interactions/PCp_CMt_intra_plts")

# PC to vSMCs

pc_vsmcs <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                               view_sel = "intra",
                               cell_pair = c("vSMCs", "PC"))

pdf("./results/sample_comparison/spatial/PC_interactions/vSMCspred_PCtar_intra_box.pdf", height = 3.5, width = 3.5)

plot(pc_vsmcs$bplot[[2]])

write_csv(pc_vsmcs$data[[2]], 
          file = "./results/sample_comparison/spatial/PC_interactions/vSMCspred_PCtar_intra_box.csv")

dev.off()


pc_vsmc_plt_dat <- bind_rows(show_useful_slides(imp_res = pc_vsmcs,
                                              pname = "vSMCs",
                                              tname = "PC",
                                              pg_name = "myogenic-enriched") %>%
                               head(),
                           show_useful_slides(imp_res = pc_vsmcs,
                                              pname = "vSMCs",
                                              tname = "PC",
                                              pg_name = "fibrotic-enriched") %>%
                             tail())


plot_visium_panels(plot_data = pc_vsmc_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/PC_interactions/vSMCsp_PCt_intra_plts")

# Endo to vSMCs

endo_PC <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                                  view_sel = "intra",
                                  cell_pair = c("Endo", "PC"))

pdf("./results/sample_comparison/spatial/PC_interactions/Endopred_PCtar_intra_box.pdf", height = 3.5, width = 3.5)

plot(endo_PC$bplot[[1]])

write_csv(endo_PC$data[[1]], 
          file = "./results/sample_comparison/spatial/PC_interactions/Endopred_PCtar_intra_box.csv")

dev.off()


endo_PC_plt_dat <- bind_rows(show_useful_slides(imp_res = endo_PC,
                                                pname = "Endo",
                                                tname = "PC",
                                                pg_name = "ischemic-enriched") %>%
                               head(),
                             show_useful_slides(imp_res = endo_PC,
                                                pname = "Endo",
                                                tname = "PC",
                                                pg_name = "fibrotic-enriched") %>%
                               tail())


plot_visium_panels(plot_data = endo_PC_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/PC_interactions/Endop_PCt_intra_plts")



# vSMCs to Fib

vsmcs_fib <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                      view_sel = "juxta_5",
                      cell_pair = c("Fib", "vSMCs"))


pdf("./results/sample_comparison/spatial/Fib_interactions/vSMCspred_Fibtar_juxta_box.pdf", height = 4.5, width = 2)

plot(vsmcs_fib$bplot[[2]])

write_csv(vsmcs_fib$data[[2]], 
          file = "./results/sample_comparison/spatial/Fib_interactions/vSMCspred_Fibtar_juxta_box.csv")

dev.off()


vsmcs_fib_plt_dat <- bind_rows(show_useful_slides(imp_res = vsmcs_fib,
                                              pname = "vSMCs",
                                              tname = "Fib",
                                              pg_name = "myogenic-enriched") %>%
                                 head(),
                           show_useful_slides(imp_res = vsmcs_fib,
                                              pname = "vSMCs",
                                              tname = "Fib",
                                              pg_name = "ischemic-enriched") %>%
                             tail())


plot_visium_panels(plot_data = vsmcs_fib_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/Fib_interactions/vSMCsp_Fibt_intra_plts")

# Immune relations

# Myeloid/Lymphoid

immune <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                                view_sel = "juxta_5",
                                cell_pair = c("Myeloid", "Lymphoid"))

pdf("./results/sample_comparison/spatial/immune_interactions/Lymphpred_Myeloidtar_intra_box.pdf", height = 4.5, width = 2)

plot(immune$bplot[[1]])

write_csv(immune$data[[1]], 
          file = "./results/sample_comparison/spatial/immune_interactions/Lymphpred_Myeloidtar_intra_box.csv")

dev.off()


immune_plt_dat <- bind_rows(show_useful_slides(imp_res = immune,
                                                  pname = "Lymphoid",
                                                  tname = "Myeloid",
                                                  pg_name = "ischemic-enriched") %>%
                                 head(),
                               show_useful_slides(imp_res = immune,
                                                  pname = "Lymphoid",
                                                  tname = "Myeloid",
                                                  pg_name = "myogenic-enriched") %>%
                                 tail())


plot_visium_panels(plot_data = immune_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/immune_interactions/Lymphp_Myeloidt_intra_plts")

# Lymphoid to CM

lymph_cm <- plot_importance_boxes(sample_importances_filt = sample_importances_filt,
                                view_sel = "intra",
                                cell_pair = c("CM", "Lymphoid"))

pdf("./results/sample_comparison/spatial/immune_interactions/Lymphpred_CMtar_intra_box.pdf", height = 3.5, width = 3.5)

plot(lymph_cm$bplot[[2]])

write_csv(lymph_cm$data[[2]], 
          file = "./results/sample_comparison/spatial/immune_interactions/Lymphpred_CMtar_intra_box.csv")

dev.off()

lymph_cm_plt_dat <- bind_rows(show_useful_slides(imp_res = lymph_cm,
                                               pname = "Lymphoid",
                                               tname = "CM",
                                               pg_name = "myogenic-enriched") %>%
                              head(),
                            show_useful_slides(imp_res = lymph_cm,
                                               pname = "Lymphoid",
                                               tname = "CM",
                                               pg_name = "ischemic-enriched") %>%
                              tail())

plot_visium_panels(plot_data = lymph_cm_plt_dat,
                   assay = "c2l",
                   out_folder = "./results/sample_comparison/spatial/immune_interactions/Lymphp_CMt_intra_plts")

