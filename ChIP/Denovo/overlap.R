setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Denovo")
denovo <- read.table("Denovo_peaks.txt", header = T, sep = "\t")
rownames(denovo) <- paste0(denovo$Chr, "_", denovo$start, "_", denovo$end)


library(GenomicRanges)

denovo <- makeGRangesFromDataFrame(denovo, keep.extra.columns = T)

# Load ATAC-seq peaks
# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

library(stringr)
atac_pks <- rownames(dea_results)
atac_peaks <- as.data.frame(str_split_fixed(atac_pks, "_", 3))
colnames(atac_peaks) <- c("Chr", "start", "end")
rownames(atac_peaks) <- atac_pks

atac_peaks <- makeGRangesFromDataFrame(atac_peaks, keep.extra.columns = T)



#----------------------------
# FIND OVERLAPPING PEAKS
#----------------------------

overlapping_peaks <- findOverlaps(atac_peaks,denovo, minoverlap = 6) # 1693 ranges 
overlap_ranges <- pintersect(atac_peaks[queryHits(overlapping_peaks)], denovo[subjectHits(overlapping_peaks)]) # shared overalping region

# Compute % occupancy per ATAC peak
atac_sizes <- width(atac_peaks[queryHits(overlapping_peaks)])
chip_occupancy <- width(overlap_ranges) / atac_sizes * 100


# save overlapping pairs # 78265 ranges
overlapping_pairs <- atac_peaks[queryHits(overlapping_peaks),]
overlapping_pairs$chip_peak <- names(denovo[subjectHits(overlapping_peaks),])
overlapping_pairs$atac_peak <- names(atac_peaks[queryHits(overlapping_peaks),])
overlapping_pairs$occupancy <- chip_occupancy

 # TRUE
head(overlapping_pairs)
summary(duplicated(overlapping_pairs$chip_peak)) 
# FALSE    TRUE 
#  1528     165

summary(duplicated(overlapping_pairs$atac_peak)) 
#FALSE    TRUE 
# 1689       4 

# Plot histogram of the occupancy
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Denovo")
pdf("occupancy_hist.pdf", width = 5, height = 5)
hist(overlapping_pairs$occupancy, main = "Occupancy distribution", xlab = "Occupancy (%)", col = "skyblue", border = "black")
dev.off()

summary(overlapping_pairs$occupancy)

# Check distribution of these peaks between differential results 
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")

# Plot sankey plot using only denovo peaks 
genes_toPlot <- overlapping_pairs$atac_peak[unique(overlapping_pairs$atac_peak)%in%unique(c(MGUS_HC_filter, SMM_HC_filter, MM_HC_filter))]

toPlot <- dea_results[genes_toPlot, c("sig_MGUS_HC","sig_SMM_HC")]
toPlot <- cbind(toPlot, dea_res_noscore[genes_toPlot, c("sig_MM_HC")])

#### PREP THE DATA ####
colnames(toPlot) <- c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC")

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Gene <- rownames(toPlot)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial
genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(MGUS_vs_HC, SMM_vs_HC, MM_vs_HC),
  names_to = "Condition",
  values_to = "Category"
)

# Order conditions
genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC"))

# Plot
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Denovo")
png("atac_dea_sanky.png", width=1500, height=1000, res=300)
ggplot(genes_long, aes(
  x = Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
dev.off()



# Check overlap with active regulatory regions | ATAC + ChIP
# Load peaks filterd by H3K27ac
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
atac_chip <- readRDS("peaks_association.RDS")

overlapping_pairs_acetyl <- overlapping_pairs[overlapping_pairs$atac_peak%in%atac_chip$peak_ID,]
summary(duplicated(overlapping_pairs_acetyl$chip_peak)) # 1521

# Check distribution of these peaks between differential results 
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")

# Plot sankey plot using only denovo peaks 
genes_toPlot <- overlapping_pairs_acetyl$atac_peak[unique(overlapping_pairs_acetyl$atac_peak)%in%unique(c(MGUS_HC, SMM_HC, MM_HC))]

toPlot <- dea_results[genes_toPlot, c("sig_MGUS_HC","sig_SMM_HC")]
toPlot <- cbind(toPlot, dea_res_noscore[genes_toPlot, c("sig_MM_HC")])

#### PREP THE DATA ####
colnames(toPlot) <- c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC")

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Gene <- rownames(toPlot)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial
genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(MGUS_vs_HC, SMM_vs_HC, MM_vs_HC),
  names_to = "Condition",
  values_to = "Category"
)

# Order conditions
genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC"))

# Plot
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Denovo")
png("atac_dea_sanky_H3K27ac.png", width=1500, height=1000, res=300)
ggplot(genes_long, aes(
  x = Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
dev.off()




# Check overlap with differentially accessible peaks 
length(MGUS_HC[MGUS_HC%in%unique(overlapping_pairs$atac_peak)]) # 487
length(SMM_HC[SMM_HC%in%unique(overlapping_pairs$atac_peak)]) # 621
length(MM_HC[MM_HC%in%unique(overlapping_pairs$atac_peak)]) # 743

summary(unique(overlapping_pairs$atac_peak)%in%unique(c(MGUS_HC, SMM_HC, MM_HC))) # 757
