#```{r}
library(cicero)
library(Seurat)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(Signac)
library(yaml)
library(Rcpp)
library(stringr)
library(WriteXLS)


AllOptions <- function(){
    parser <- OptionParser()
    parser <- add_option(parser, c("-s", "--SAMPLE"), type="character", default="CK166",
                    help="sample name for ATAC integration [default %default]",
                    metavar="character")
    return(parser)
}
parser <- AllOptions()
pa <- parse_args(parser)

atac_name   = pa$SAMPLE
atac_to_rna = setNames( c("CK158", "CK159", "CK160",  "CK162", "CK163", "CK164","CK165"), 
                    c("CK166", "CK167", "CK168",  "CK170", "CK171", "CK173", "CK174"))
rna_name = atac_to_rna[atac_name]


atac_dir <- paste0("../ATAC_Single_Sample_V2/data/", atac_name)
###--------------unionPeaks--------------------------

u_atac_name <- paste0(atac_name, "_unionPeaks.Rds")


set.seed(1)

getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
      #seqnames=#paste0("chr",seqnames(geneGTF1)),
      seqnames=seqnames(geneGTF1),
      ranges=ranges(geneGTF1),
      strand=strand(geneGTF1),
      gene_name=geneGTF1$gene_name,
      gene_id=geneGTF1$gene_id
    #) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
    ) %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}


fixA <- "center"
fixB <- "start"
associationWindow <- 2 * 250*10^3 + 1 #+-250 Kb
corCutOff <- 0.35 #Pearson Correlation Cutoff
fdrCutOff <- 0.1 #FDR Cutoff
distCutOff <- 2500 #Min Dist to TSS


#load(file="save/rna_a_ZKY095.Rdata")
#load(file="save/atac_a_ZKY095.Rdata")

rna_dir <- "../scRNA_filtered/" 

rna_a <- readRDS(file=file.path(rna_dir, sprintf("%s.filtered.annotated.rds", rna_name)))

gtf <- getGeneGTF("static/genes.gtf")

DefaultAssay(rna_a) <- "RNA"

inter_gene_name <- intersect(rownames(rna_a), gtf@elementMetadata$gene_name)
rna_a <- subset(rna_a, features=inter_gene_name)


seRNA <- SummarizedExperiment(
  assays = SimpleList(counts=rna_a@assays$RNA@counts)
)

rm(rna_a)
gc()


gtfMatch <- gtf[na.omit(match(rownames(seRNA), gtf$gene_name))]
names(gtfMatch) <- rownames(seRNA)
rowRanges(seRNA) <- gtfMatch



atac <- readRDS(file=file.path(atac_dir, u_atac_name))
DefaultAssay(atac) <- "peaks"

atac_mtx <- GetAssayData(atac, slot = "counts", assay = "peaks")

#DefaultAssay(atac_a) <- "Peaks"

peaks <- rownames(atac_mtx)
peak.ranges <- StringToGRanges(regions =peaks , sep = c(":", "-"))
seATAC <- SummarizedExperiment(
                        assays <- SimpleList(atac_mtx),
                        colData = atac@meta.data,
                        rowRanges = peak.ranges)


rm(atac_mtx)
gc()
#assay(seATAC) <- log2(edgeR::cpm(assay(seATAC))/100+1)
#assay(seRNA) <- log2(edgeR::cpm(assay(seRNA))/100+1)


seRNAWindow <- resize(rowRanges(seRNA), width = 1, fix = fixB) %>%
  {suppressWarnings(resize(., width = associationWindow, fix = "center"))} %>% trim(.)

#Keep only seATAC within association window
seATAC <- seATAC[unique(queryHits(findOverlaps(resize(seATAC,1,fixA), seRNAWindow, ignore.strand = TRUE)))]

#Getting distances
message("Getting Distances...")
o <- findOverlaps(seRNAWindow, resize(rowRanges(seATAC),1,fixA), ignore.strand = TRUE)


#Get Distance from Fixed point A B correct for minus stranded
mcols(o)$distance <- start(resize(rowRanges(seATAC),1,fixA))[subjectHits(o)] - start(resize(rowRanges(seRNA),1,fixB))[queryHits(o)]
mcols(o)$distance[which(as.character(strand(rowRanges(seRNA)))[queryHits(o)]=="-")] <- -1*mcols(o)$distance[which(as.character(strand(rowRanges(seRNA)))[queryHits(o)]=="-")]

#Add other info
o <- DataFrame(o)
colnames(o) <- c("B","A","distance")
o <- o[,c("A","B","distance")]

#```




#```{r}
#Get GTF
#gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
strand(tssRNA) <- "*" 
peakLinks <- rowRanges(seATAC)[o[,1]]
geneLinks <- rowRanges(seRNA) %>% resize(1, "start") %>% {.[o[,2]]}

mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")

rm(seRNA, seATAC)
gc()


#peakName geneName peak_chr peak_start peak_end gene_chr gene_start gene_strand gene_name gene_id exonLength distance
out_list <- peakLinks@elementMetadata@listData
out_list$geneName <- out_list$gene_name
out_list$peak_chr <- sapply(str_split(out_list$peakName, "_"), function(x) x[1]) 
out_list$peak_start <- sapply(str_split(out_list$peakName, "_"), function(x) x[2])
out_list$peak_end <- sapply(str_split(out_list$peakName, "_"), function(x) x[3])
out_list$peakName <- sub("_", ":",  out_list$peakName)
out_list$peakName <- sub("_", "-",  out_list$peakName)


out_list <- out_list[c("peakName", "geneName", "peak_chr",
                          "peak_start", "peak_end", "gene_chr",
                          "gene_start", "gene_strand", "gene_name",
                          "gene_id", "exonLength", "distance")]

#peak    gene    distance    peak_type
tsv_list <- out_list[c("peakName", "geneName", "distance")]
tsv_df <- as.data.frame(tsv_list)
tsv_df <- rename(tsv_df, c("peakName" = "peak",
                           "geneName" = "gene"))
write.table(tsv_df, file=sprintf("save/peak2gene_putative-simple_%s.tsv", atac_name), sep="\t", quote=F, row.names = F)

write.table(as.data.frame(out_list), file=sprintf("save/peak2gene_putative_%s.tsv", atac_name), sep="\t", quote=F, row.names = F)



#WriteXLS(as.data.frame(out_list),
#         file.path("data", "peak2gene_putative.xlsx"),
#         SheetNames = "peak2gene_putative")



#```
