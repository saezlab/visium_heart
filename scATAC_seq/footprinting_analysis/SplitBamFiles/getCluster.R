library(stringr)
library(Seurat)
library(plyr)


sample_list <- c("CK166", "CK167", "CK168", "CK170", "CK171", "CK173", "CK174")
for (sample in sample_list) {
    obj_filename <- sprintf("../../ATAC_Single_Sample_V2/data/%s/%s_unionPeaks.Rds", 
                            sample, sample)
    obj <- readRDS(obj_filename)
    df <- obj@meta.data
    df$Barcode <- rownames(df)
    df <- subset(df, select = c("Barcode", "peaks_snn_res.0.5"))
    
    colnames(df) <- c("Barcode", "Cluster")
    write.table(df, file = paste0("./Clusters/", sample, ".txt"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
}
