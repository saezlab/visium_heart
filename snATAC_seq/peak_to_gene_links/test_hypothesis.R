library(dplyr)
library(doParallel)
registerDoParallel(cores=8)


## please add adjust pval to it


atac_names = c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174")
#atac_names = c("CK167", "CK168")#, "CK171", "CK173", "CK174")
#atac_names = c("CK166")
#atac_names = c("CK166")#, "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174")

foreach(atac_name = atac_names) %dopar%{
    message("processing ", atac_name, date())
    df <- read.csv(file=sprintf("save/corr_%s.tsv", atac_name), sep="\t", stringsAsFactors = F)
    fn = sprintf("save/%s_null_hypothesis.Rds", atac_name)
    null_dist_list <- readRDS(file=fn)
    stt_vars <- lapply(null_dist_list, function(lst) return(c("mean"=mean(lst), "sd"=sd(lst))))
    
    #mean_null_coor = 0.1
    #sd_null_corr = 0.2
    df$pval <- sapply(1:nrow(df), function(x){
              chromosome <- strsplit(df[x,]$peak, ":")[[1]][1]
              mean_null_corr <- stt_vars[[chromosome]]["mean"]
              sd_null_corr <- stt_vars[[chromosome]]["sd"]
              2*pnorm(-abs(((df[x,]$corr - mean_null_corr) / sd_null_corr)))
      })
    
    df$padj <- p.adjust(df$pval, method="BH")
    
    write.table(df, file=sprintf("save/corr_pval_%s.tsv",atac_name), sep = "\t", quote = F, row.names = F)

    message("finished", atac_name, date())
}

