library(dplyr)
library(VennDiagram)
#library(RColorBrewer)
library(tidyverse)
library(optparse)
#futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


parser <- OptionParser()
parser <- add_option(parser, c("-f", "--filter"), type="character", default="0.01",
                    help="filter pval [default %default]",
                    metavar="character")

parser <- add_option(parser, c("-s", "--sample"), type="character", default="CK166",
                    help="atac name [default %default]",
                    metavar="character")

pa = parse_args(parser)

filter_pval <- as.numeric(pa$filter)
atac_name <- pa$sample

df <- read.csv(file=sprintf("save/corr_pval_%s.tsv", atac_name), sep = "\t")

fdf <- df %>% filter(corr>0)    
fdf$Sig <- ifelse(fdf$padj<0.05, "TRUE", "FALSE")

fdf <- fdf %>% filter(fdf$Sig == "TRUE")
write.table(fdf, file=sprintf("save/corr_pval_final_%s.tsv", atac_name),
                    sep = "\t", quote = F, row.names = F)

fdf <- fdf %>% filter(fdf$padj<filter_pval)
fdf <- fdf %>% filter(abs(distance) > 2000)
write.table(fdf, file=sprintf("save/corr_pval_final_%s_%.2f.tsv", atac_name, filter_pval),
                    sep = "\t", quote = F, row.names = F)


#df <- read.csv(file=sprintf("save/corr_pval_final_%s_0.01.tsv", atac_name), sep = "\t")
#df <- df %>% filter(df$corr > 0.5)
write.table(df %>% filter(df$corr > 0.3), 
                file=sprintf("save/corr_pval_final_%s_%.2f_0.3.tsv",atac_name, filter_pval), 
                sep = "\t", quote = F, row.names = F)

write.table(df %>% filter(df$corr > 0.5), 
                file=sprintf("save/corr_pval_final_%s_%.2f_0.5.tsv",atac_name, filter_pval), 
                sep = "\t", quote = F, row.names = F)


