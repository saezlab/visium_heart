get_celltypes <- function(plots_opt, atac_name){

if(plots_opt == ""){

    cluster_anno <- switch(atac_name,
"CK166" = c("8" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "2" = "Endothelial cells 1",
  "6" = "Endothelial cells 2 POSTN+",
  "3" = "Fibroblasts 1 COL15A1+",
  "7" = "Macrophages 1 CD163+",
  "5" = "Pericytes",
  "4" = "Pericytes"),
  
 "CK167" = c("5" = "Cardiomyocytes",
  "3" = "Endothelial cells 1",
  "8" = "Endothelial cells 1",
  "2" = "Endothelial cells 2 POSTN+",
  "1" = "Fibroblasts",
  "4" = "Fibroblasts",
  "7" = "Macrophages",
  "6" = "Pericytes"),


"CK168" = c("8" = "Cardiomyocytes 1",
  "6" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "3" = "Macrophages",
  "4" = "Macrophages",
  "7" = "Pericytes EGFLAM+"),

"CK170" = c("5" = "Cardiomyocytes",
  "2"= "Cardiomyocytes",
  "4" = "Damaged endothelial cells",
  "1" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "6" = "Fibroblasts"),

"CK171" = c("1" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "6" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "5" = "Fibroblasts",
  "4" = "Macrophages",
  "9" = "Neuronal cells",
  "2" = "Pericytes",
  "3" = "Vascular smooth muscle cells"),

"CK173" = c("3" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "6" = "Endothelial cells 1",
  "5" = "Endothelial cells 1",
  "8" = "Fibroblasts 1",
  "9" = "Fibroblasts 1",
  "7" = "Fibroblasts 2 SCARA5+",
  "10" = "Macrophages",
  "4" = "Pericytes"),

"CK174" = c("2" = "Cardiomyocytes 1",
  "6" = "Cardiomyocytes 2",
  "3" = "Endothelial cells",
  "4" = "Fibroblasts",
  "5" = "Macrophages",
  "1" = "Pericytes"))
}

if(plots_opt == "_0.01"){
    # 4,6,5,9,3,1,2,8,7

   cluster_anno <- switch(atac_name,
   "CK166" = c("1" = "Cardiomyocytes",
  "5" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "8" = "Endothelial cells 2 POSTN+",
  "2" = "Fibroblasts 1 COL15A1+",
  "7" = "Macrophages 1 CD163+",
  "3" = "Pericytes",
  "6" = "Pericytes"),
  
 "CK167" = c("4" = "Cardiomyocytes",
  "8" = "Endothelial cells 1",
  "7" = "Endothelial cells 2 POSTN+",
  "2" = "Endothelial cells 3 PLVAP+",
  "1" = "Fibroblasts",
  "3" = "Fibroblasts",
  "6" = "Macrophages",
  "5" = "Pericytes"),

"CK168" = c("8" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 1",
  "3" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "4" = "Fibroblasts 1",
  "6" = "Macrophages",
  "7" = "Pericytes EGFLAM+"),

"CK170" = c("3" = "Cardiomyocytes",
  "1" = "Cardiomyocytes",
  "2" = "Damaged endothelial cells",
  "5" = "Endothelial cells 1",
  "4" = "Fibroblasts",
  "6" = "Fibroblasts"),

"CK171" = c("1" = "Cardiomyocytes",
  "7" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "5" = "Macrophages",
  "6" = "Neuronal cells",
  "2" = "Pericytes",
  "8" = "Pericytes",
  "9" = "Vascular smooth muscle cells"),

"CK173" = c("1" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "9" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "6" = "Fibroblasts 1",
  "5" = "Fibroblasts 2 SCARA5+",
  "10" = "Macrophages",
  "3" = "Pericytes"),

"CK174" = c("6" = "Cardiomyocytes 1",
  "2" = "Endothelial cells",
  "3" = "Fibroblasts",
  "5" = "Fibroblasts",
  "4" = "Macrophages",
  "1" = "Pericytes")
        )

    
}

if(plots_opt == "_0.01_0.3"){
    #3,1,4,6,8,2,7,9,5

cluster_anno <- switch(atac_name,
"CK166" = c("4" = "Cardiomyocytes",
  "3" = "Endothelial cells 1",
  "8" = "Endothelial cells 1",
  "7" = "Endothelial cells 2 POSTN+",
  "1" = "Fibroblasts 1 COL15A1+",
  "2" = "Fibroblasts 1 COL15A1+",
  "5" = "Macrophages 1 CD163+",
  "6" = "Pericytes"),
  
  
  
"CK167" = c("5" = "Cardiomyocytes",
  "1" = "Endothelial cells 1",
  "7" = "Endothelial cells 1",
  "3" = "Fibroblasts",
  "2" = "Fibroblasts",
  "6" = "Macrophages",
  "8" = "Neuronal cells",
  "4" = "Pericytes"),


"CK168" = c("8" = "Cardiomyocytes 1",
  "5" = "Cardiomyocytes 1",
  "3" = "Cardiomyocytes 2",
  "1" = "Endothelial cells 1",
  "2" = "Fibroblasts 1",
  "6" = "Macrophages",
  "4" = "Pericytes EGFLAM+",
  "7" = "Pericytes EGFLAM+"),


"CK170" = c("1" = "Cardiomyocytes",
  "3" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "6" = "Endothelial cells 1",
  "2" = "Fibroblasts",
  "5" = "Fibroblasts"),



"CK171" = c("2" = "Cardiomyocytes",
  "5" = "Endothelial cells 1",
  "9" = "Endothelial cells 1",
  "1" = "Fibroblasts",
  "6" = "Macrophages",
  "7" = "Neuronal cells",
  "4" = "Pericytes",
  "3" = "Pericytes",
  "8" = "Vascular smooth muscle cells"),



"CK173" = c("1" = "Cardiomyocytes",
  "7" = "Cardiomyocytes",
  "8" = "Cardiomyocytes",
  "2" = "Cardiomyocytes",
  "4" = "Endothelial cells 1",
  "6" = "Endothelial cells 1",
  "5" = "Fibroblasts 1",
  "10" = "Fibroblasts 1",
  "9" = "Macrophages",
  "3" = "Pericytes"),


"CK174" = c("2" = "Cardiomyocytes 1",
  "4" = "Cardiomyocytes 2",
  "3" = "Endothelial cells",
  "6" = "Fibroblasts",
  "5" = "Macrophages",
  "1" = "Pericytes")

)
}

	return(cluster_anno)
}

suppressPackageStartupMessages(library(optparse))


parser <- OptionParser()



parser <- add_option(parser, c("-p", "--plots"), type="character", default="",
                    help="plots for different pval_option [default %default]",
                    metavar="character")

pa = parse_args(parser)

plots_opt <-  ifelse(pa$plots == "",  "",  paste0("_", pa$plots))
print(plots_opt)
dir.create(sprintf("plots%s", plots_opt))

sample_vec <- c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174")  


pb = txtProgressBar(min = 0, max = length(sample_vec), initial = 0) 
for(a_sample in sample_vec){
	a_name = sprintf("save/detail_corr_pval_final_%s%s.tsv" , a_sample, plots_opt)
	print(a_name)
	df <- read.csv(file=a_name , sep="\t")
    df[, paste0(a_sample, "_clusters")] <- NA

	celltype.list <- get_celltypes(plots_opt, a_sample)
	number.clusters <- length(celltype.list)

    fname <- sprintf("plots%s/merge_info_%s_%d.Rds", plots_opt, a_sample, number.clusters)    
    c_list <- readRDS(fname)
    setTxtProgressBar(pb, which(sample_vec == a_sample) )
    message(a_sample, "\t",date())
    for(nm in names(c_list)){
        a_list <- c_list[[nm]]
        nums <- length(a_list[[1]])
        message(nm, date())
        for(i in 1:nums){
           peak <- a_list[[1]][i]
           gene <- a_list[[2]][i]
		   celltype = celltype.list[as.character(nm)]
           df[df$peakName==peak & df$geneName == gene,][, paste0(a_sample, "_clusters")]<-celltype
        }
    }
    fname <- sprintf("save/detail_corr_pval_final_%s%s_with_annotation.tsv", a_sample,plots_opt)
    write.table(df, file = fname, sep="\t", row.names=F,quote=FALSE)
}



