###########################################
####         ATAC-Seq pipeline         ####
####          Peak annotation          ####
####              HUMAN                ####
###########################################

# Library required
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)


# Load consensus peaks
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")


## Annotation
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene


#>>>>>>>>>>>>>>>>> Create the GRanges of consensus peaks

### Filtered peaks
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA/")
load("voom_to_dea.RData")
peaks <- rownames(norm_data)
consensus_peaks <- consensus_peaks[peaks,] 

# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")
counts <- readRDS("ATAC_raw_counts.RDS") # Load the count data
dim(counts) # 51685 x 49


# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == rownames(consensus_peaks))
# TRUE


#compute the mean coverage
count_table <- counts
mean_cov <- rowMeans(count_table)

# Read consensus peaks
peaks <- consensus_peaks
peaks$coverage <- mean_cov

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)


##Coverage plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Peaks")
png("peaks_over_chr.png")
covplot(peak_granges, title = "ATAC-Seq Peaks over Chromosomes", weightCol="coverage")
dev.off()


# Profile of ATAC peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak_granges, windows=promoter)

# Heatmap of ATAC binding to TSS regions
png("binding_TSS.png")
tagHeatmap(tagMatrix)
dev.off()


# Average Profile of ATAC peaks binding to TSS region
png("average_binding_TSS.png")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()


# Peak Annotation
peakAnno <- annotatePeak(peak_granges, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db", level="gene")

peak_annotation <- as.data.frame(peakAnno)

peak_annotation$annotation2 <- peak_annotation$annotation
peak_annotation$annotation2[grep("Promoter",peak_annotation$annotation)] <- "Promoter"
peak_annotation$annotation2[grep("Intron",peak_annotation$annotation)] <- "Intron"
peak_annotation$annotation2[grep("Distal Intergenic",peak_annotation$annotation)] <- "Distal Intergenic"
peak_annotation$annotation2[grep("Exon",peak_annotation$annotation)] <- "Exon"

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Peaks")
saveRDS(peak_annotation, file="peak_annotation.RDS") 


