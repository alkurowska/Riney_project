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
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")


## Annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#>>>>>>>>>>>>>>>>> Create the GRanges of consensus peaks

### Filtered peaks
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
load("voom_to_dea.RData")
peaks <- rownames(norm_data)


# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata") 
counts <- as.data.frame(Counts$counts) # 157817 x  248

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == consensus_peaks$GeneID)
# TRUE
rownames(counts) <- consensus_peaks$name # change rownames of the count data

# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
counts <- counts[,-which(colnames(counts) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]

# change colnames of the count data
# Remove everything after "_S"
colnames(counts) <- gsub("_S.*", "", colnames(counts))


# samples to keep
dim(counts)
counts <- counts[,colnames(counts) %in% colnames(norm_data)] # 157817 x  185
counts <- counts[,colnames(norm_data)] # 185

consensus_peaks <- consensus_peaks[peaks,]
counts <- counts[consensus_peaks$name,] # 142136 x  185
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
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
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
peakAnno <- annotatePeak(peak_granges, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db", level="gene")

peak_annotation <- as.data.frame(peakAnno)

peak_annotation$annotation2 <- peak_annotation$annotation
peak_annotation$annotation2[grep("Promoter",peak_annotation$annotation)] <- "Promoter"
peak_annotation$annotation2[grep("Intron",peak_annotation$annotation)] <- "Intron"
peak_annotation$annotation2[grep("Distal Intergenic",peak_annotation$annotation)] <- "Distal Intergenic"
peak_annotation$annotation2[grep("Exon",peak_annotation$annotation)] <- "Exon"

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
saveRDS(peak_annotation, file="peak_annotation.RDS") 


