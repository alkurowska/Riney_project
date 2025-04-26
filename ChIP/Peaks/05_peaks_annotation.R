###########################################
####         ATAC-Seq pipeline         ####
####          Peak annotation          ####
####              HUMAN                ####
###########################################

# Library required
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load consensus peaks
setwd("/ibex/user/kurowsaa/RINEY/NEW/ChIP/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

## Annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#>>>>>>>>>>>>>>>>> Create the GRanges of consensus peaks

### Filtered peaks
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/")
load("norm_to_voom.RData")
peaks <- rownames(y$counts)


# LOAD THE RAW DATA

# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix")
chip_data <- readRDS("CHIP_Info.RDS")
counts <- as.data.frame(chip_data$counts) # 94852 x  136

# Consensus peaks
consensus_peaks <- chip_data$genes
consensus_peaks$name <- paste0(consensus_peaks$Chr,"_",consensus_peaks$Start,"_",consensus_peaks$End)
rownames(consensus_peaks) <- consensus_peaks$name

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == consensus_peaks$GeneID)
# TRUE
rownames(counts) <- consensus_peaks$name # change rownames of the count data


# samples/peaks to keep
dim(counts)
counts <- counts[rownames(counts)%in%peaks,colnames(counts) %in% colnames(y$counts)] # 93674 x 128

consensus_peaks <- consensus_peaks[peaks,]

table(rownames(counts) == consensus_peaks$name)
rownames(counts) <- consensus_peaks$name

#compute the mean coverage
count_table <- counts
mean_cov <- rowMeans(count_table)

# Read consensus peaks
peaks <- consensus_peaks
peaks$coverage <- mean_cov

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)


##Coverage plot
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Peaks")
png("peaks_over_chr.png")
covplot(peak_granges, title = "ChIP-Seq Peaks over Chromosomes", weightCol="coverage")
dev.off()


# Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak_granges, windows=promoter)

# Heatmap of ChIP binding to TSS regions
png("binding_TSS.png")
tagHeatmap(tagMatrix)
dev.off()


# Average Profile of ChIP peaks binding to TSS region
png("average_binding_TSS.png")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()


# Peak Annotation
peakAnno <- annotatePeak(peak_granges, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db", level="gene")

peak_annotation <- as.data.frame(peakAnno)

peak_annotation$annotation2 <- peak_annotation$annotation
peak_annotation$annotation2[grep("Promoter",peak_annotation$annotation)] <- "Promoter"
peak_annotation$annotation2[grep("Intron",peak_annotation$annotation)] <- "Intron"
peak_annotation$annotation2[grep("Distal Intergenic",peak_annotation$annotation)] <- "Distal Intergenic"
peak_annotation$annotation2[grep("Exon",peak_annotation$annotation)] <- "Exon"

# Save data
saveRDS(peak_annotation, file="peak_annotation.RDS") 
