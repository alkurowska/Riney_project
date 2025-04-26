#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: chromVAR_TF.R          #
#***********************************#


# https://greenleaflab.github.io/chromVAR/articles/Introduction.html
# https://bioconductor.org/packages/3.3/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment

# loading required packages
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(purrr)
library(TFBSTools)
library(JASPAR2024)

# First we need to construct SummarizedExperiment object 

## Getting the required files
# 1- Set of OCR
# Consensus Peaks 
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# 2- OCR counts and metadata
setwd("/ibex/user/kurowsaa/RINEY/human/HC_MGUS/ATAC")
load("Rsubread_Counts_ATAC_Study1.Rdata")
my_counts_matrix <- Counts$counts
peaks <- Counts$annotation
rownames(peaks) <- paste(peaks$Chr, peaks$Start, peaks$End, sep = "_")

# Samples to keep

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(my_counts_matrix) == consensus_peaks$GeneID)
# TRUE
rownames(my_counts_matrix) <- consensus_peaks$name # change rownames of the count data

# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]

# change colnames of the count data
# Remove everything after "_S"
colnames(my_counts_matrix) <- gsub("_S.*", "", colnames(my_counts_matrix))

# load metadata
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
load("norm_to_voom.RData")

my_counts_matrix <- my_counts_matrix[,colnames(my_counts_matrix)%in%fish_atac$`Sample ATACseq`]
dim(my_counts_matrix) # 157817  x  185

fish_atac <- fish_atac[fish_atac$`Sample ATACseq` %in% colnames(my_counts_matrix),]
dim(fish_atac) # 185 x 17

peaks <- peaks[rownames(y$counts),]
my_counts_matrix <- my_counts_matrix[rownames(peaks),]

table(rownames(my_counts_matrix) == rownames(peaks))


# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks, keep.extra.columns = T)
peak_granges

## Creating the object 

fragment_counts <- SummarizedExperiment(assays=list(counts=my_counts_matrix),
                                        rowRanges=peak_granges, colData=fish_atac)

# Getting GC content of peaks

fragment_counts2 <- addGCBias(fragment_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg38)

#Filtering peaks
counts_filtered <- filterPeaks(fragment_counts2, non_overlapping = TRUE)


#Get motifs and what peaks contain motifs
## JASPAR 2024
#https://academic.oup.com/nar/article/46/D1/D252/4616875
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks/Jaspar_motifs")

#Read the motifs
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
motifs24 <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))

# Save motifs names 
names <- as.data.frame(name(motifs24))
save(names, file="JASPAR2024.RData")


# keep only expressed
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
load("norm_to_voom.RData")
genes <- unique(gene_anno[gene_anno$gene_id%in%rownames(d1_norm$counts),]$gene_name)
TFs <- unique(names[,1][names[,1]%in%genes]) #453

motifs24 <- motifs24[which(names[,1]%in%TFs)]
# Get motif matches in peaks
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks/Jaspar_motifs")
motif_ix <- matchMotifs(motifs24, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#Save as RData
save(motif_ix, file = "Jaspar_motifs_matches.RData")  

# Get motif positions within peaks for motifs in peaks 

jaspar_ranges <- matchMotifs(motifs24, peak_granges, genome = "hg38", out = "positions") 
for(i in 1:length(motifs24)){
  names(jaspar_ranges)[i] <- unname(name(motifs24[i]))
}
save(jaspar_ranges, file = "Jaspar_motifs_positions.RData")
##>>>>>>>>>>>>>>>>>>>>> Relate each OCR to motif


jaspar_names <- names(jaspar_ranges)

ocr_jaspar_matrix <- matrix(0,nrow=nrow(consensus_peaks),ncol=length(jaspar_names))
colnames(ocr_jaspar_matrix) <- jaspar_names
names(peak_granges) <- paste0(consensus_peaks$Chr,"_", consensus_peaks$Start, "_", consensus_peaks$End)
rownames(ocr_jaspar_matrix) <- names(peak_granges)

for(i in 1:length(jaspar_names)){
  cat(i," - ")
  #select granges for target motif
  jaspar_gr_motif <- jaspar_ranges[[i]]
  #ranges overlaping
  ocr_jaspar_association <- findOverlaps(jaspar_gr_motif,peak_granges,type = "within")
  if(length(ocr_jaspar_association)>0){
    ocr_jaspar_matrix[unique(ocr_jaspar_association@to),i] <- 1    
  }  
}

colnames(ocr_jaspar_matrix) <- jaspar_names
write.table(ocr_jaspar_matrix, "ocr_jaspar_matrix.txt", sep = "\t", col.names = T, row.names = T)

