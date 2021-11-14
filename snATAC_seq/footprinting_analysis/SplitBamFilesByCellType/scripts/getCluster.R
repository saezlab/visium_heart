library(stringr)
library(ArchR)

obj.atac <- readRDS("../../../DataIntegration/data/VisiumHeart/snATAC.annotated.Rds")

df <- obj.atac@meta.data %>%
    as.data.frame()

df <- subset(df, select = c("Sample", "cell_type"))
df$Barcode <- rownames(df)

# change the barcode suffex
df$Barcode <- stringr::str_split_fixed(df$Barcode, "#", 2)[, 2]


for (sample in unique(df$Sample)) {
    df_sub <- subset(df, Sample == sample)
    df_sub$Sample <- NULL
    df_sub <- df_sub[, c("Barcode", "cell_type")]
    
    write.table(df_sub, file = glue::glue("{sample}.txt"), 
                row.names = FALSE, sep = "\t", quote = FALSE)

}
