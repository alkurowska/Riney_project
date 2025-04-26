
#Log cpm - tmm - no batch - ATAC data
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
atac_matrix <- read.table("log_cpm_noFRIP.txt", header=T, sep="\t")  
colnames(atac_matrix) <- gsub("X", "", colnames(atac_matrix))
dim(atac_matrix) # 142136    185

#Log cpm - tmm - no batch - ChIP data
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix")
chip_matrix <- read.table("log_cpm_noFRIP.txt", header=TRUE)
colnames(chip_matrix) <- gsub("X", "", colnames(chip_matrix))
dim(chip_matrix) # 93674   197

# MATCH ATAC & ChIP names
colnames(chip_matrix)[!colnames(chip_matrix) %in% colnames(atac_matrix)]

# Paired samples between ATAC and ChIP
common <- intersect(colnames(atac_matrix), colnames(chip_matrix))
atac_common <- atac_matrix[,common]
chip_common <- chip_matrix[,common]

library(stringr)

# split by _ and conver to a data frame with 3 columns
atac_peaks <- as.data.frame(str_split_fixed(rownames(atac_common), "_", 3))
colnames(atac_peaks) <- c("chr", "start", "end")
rownames(atac_peaks) <- rownames(atac_common)

library(GenomicRanges)
atac_peaks <- makeGRangesFromDataFrame(atac_peaks, keep.extra.columns=TRUE)

chip_peaks <- as.data.frame(str_split_fixed(rownames(chip_common), "_", 3))
colnames(chip_peaks) <- c("chr", "start", "end")
rownames(chip_peaks) <- rownames(chip_common)
chip_gr <- makeGRangesFromDataFrame(chip_peaks, keep.extra.columns=TRUE)


# FIND OVERLAPS with ChIP peaks
# Find overlaps
#-------------------------------------
overlaps <- findOverlaps(atac_peaks, chip_gr, minoverlap=6) # 77616
overlap_ranges <- pintersect(atac_peaks[queryHits(overlaps)], chip_gr[subjectHits(overlaps)]) # shared overalping region

# Compute % occupancy per ATAC peak
atac_sizes <- width(atac_peaks[queryHits(overlaps)])
chip_occupancy <- width(overlap_ranges) / atac_sizes * 100


# save overlapping pairs # 78265 ranges
overlapping_pairs <- atac_peaks[queryHits(overlaps),]
overlapping_pairs$chip_peak <- names(chip_gr[subjectHits(overlaps),])
overlapping_pairs$atac_peak <- names(atac_peaks[queryHits(overlaps),])
overlapping_pairs$occupancy <- chip_occupancy

 # TRUE
head(overlapping_pairs)
summary(duplicated(overlapping_pairs$chip_peak)) 
# FALSE    TRUE 
# 69140    8476

summary(duplicated(overlapping_pairs$atac_peak)) 
#FALSE    TRUE 
#77361     255 

# Plot histogram of the occupancy
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
pdf("occupancy_hist.pdf", width = 5, height = 5)
hist(overlapping_pairs$occupancy, main = "Occupancy distribution", xlab = "Occupancy (%)", col = "skyblue", border = "black")
dev.off()

summary(overlapping_pairs$occupancy)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.2077  84.1308 100.0000  88.6022 100.0000 100.0000  
#-------------------------------------

# Keep only the overlapping peaks
chip <- chip_common[rownames(chip_common) %in% overlapping_pairs$chip_peak,]
atac <- atac_common[rownames(atac_common) %in% overlapping_pairs$atac_peak,]


chip <- as.matrix(chip)
atac <- as.matrix(atac)

dim(chip) # 69140   128
dim(atac) # 77361   128


pairs <- cbind(overlapping_pairs$chip_peak, overlapping_pairs$atac_peak)
pairs <- as.data.frame(pairs)
colnames(pairs) <- c("chip_peak", "atac_peak")

##>>>>>>>>>>>>>>>>>>>>> Correlation
peak_association <- list()

for(i in 1:nrow(pairs)){ #for each atac peak i
# select a chip peak i in a pair
  chip_peak <- pairs[i,]$chip_peak
  atac_peak <- pairs[i,]$atac_peak

  cat("peak ",i,"-peaks ","\n")
  cor_p <- vector() #p-value
  cor_r <- vector() #correlation
  chip.id <- vector() #chip ID
  peak.id <- vector() #peak ID
  cc <- cor.test(chip[chip_peak,],atac[atac_peak,],method="spearman")
  cor_p <- c(cor_p,cc$p.value) 
  cor_r <- c(cor_r,cc$estimate)
  chip.id <- c(chip.id, chip_peak)
  peak.id <- c(peak.id, atac_peak)

  cor_output <- data.frame(peak=atac_peak,spearman_pval=cor_p,spearman_rho=cor_r, chip_ID=chip.id, peak_ID=peak.id)

  peak_association[[i]] <- cor_output
}


##Let's see the signifcance of the peaks 
data <- do.call("rbind", peak_association)
data$spearman_adj.pval <- p.adjust(data$spearman_pval,method="fdr")

data$occupancy <- overlapping_pairs$occupancy

#Save results
write.table(data, "peaks_association.txt", sep="\t", dec=".", quote=FALSE, row.names=FALSE, col.names = TRUE)
saveRDS(data, file = "peaks_association.RDS")

summary(data$spearman_rho)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2911  0.1997  0.3074  0.3159  0.4249  0.8882 

dim(data[data$spearman_adj.pval<0.05,]) # 60882
dim(data[data$spearman_adj.pval<0.01,]) # 51464

# keep only significant peaks
data_sig <- data[data$spearman_adj.pval<0.05,]

summary(data_sig$spearman_rho)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.2911  0.2709  0.3536  0.3751  0.4581  0.8882 




toPlot <- cbind(data$occupancy, data$spearman_rho)
toPlot <- as.data.frame(toPlot)
colnames(toPlot) <- c("occupancy", "spearman_rho")
toPlot$significant <- 0
toPlot$significant[data$spearman_adj.pval < 0.05] <- 1

library(ggplot2)
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho, color = factor(significant))) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  # change legend labels
  scale_color_manual(values = c("black", "skyblue"),
                     labels = c("Not significant", "Significant")) +
  guides(color = guide_legend(title = "p.adj < 0.05")) + 
  # legend title
  theme_minimal()
dev.off()


# Density plot of rho

toPlot <- data$spearman_rho
library(reshape2)
toPlot_melted <- melt(toPlot)

# Create the density plot

p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot.png", p, width = 6, height = 4, dpi = 300)

# Density plot of rho significant
# Get the significant

sig <- data[data$spearman_adj.pval < 0.05,]

toPlot <- sig$spearman_rho

toPlot_melted <- melt(toPlot)

# Create the density plot

p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho of significant correlations",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot_sig.png", p, width = 6, height = 4, dpi = 300)


length(unique(data_sig$peak_ID)) # 60709 unique atac peaks
length(unique(data_sig$chip_ID)) # 55165 unique chip peaks


# Check distribution of these peaks between differential results 
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")

MM_HC <- MM_HC_filter
SMM_HC <- SMM_HC_filter
MGUS_HC <- MGUS_HC_filter

summary(data_sig$peak_ID%in%MM_HC) # 13790
summary(data_sig$peak_ID%in%SMM_HC) # 8264
summary(data_sig$peak_ID%in%MGUS_HC) # 5039


toPlot <- data_sig[data_sig$peak_ID%in%MM_HC,]
toPlot <- as.data.frame(toPlot)

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_MM.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho)) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  # legend title
  theme_minimal()
dev.off()



toPlot <- data_sig[data_sig$peak_ID%in%SMM_HC,]
toPlot <- as.data.frame(toPlot)

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_SMM.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho)) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  # legend title
  theme_minimal()
dev.off()



toPlot <- data_sig[data_sig$peak_ID%in%MGUS_HC,]
toPlot <- as.data.frame(toPlot)

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_MGUS.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho)) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  # legend title
  theme_minimal()
dev.off()


# How many up and down?

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")


summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])) # 7342 UP
summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])) # 6448 DOWN

summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_SMM_HC==1,])) # 5791 UP
summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_SMM_HC==-1,])) # 4440 DOWN

summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MGUS_HC==1,])) # 3185 UP
summary(data_sig$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MGUS_HC==-1,])) # 2263 DOWN



toPlot <- data_sig[data_sig$peak_ID%in%MM_HC,]
toPlot <- as.data.frame(toPlot)
toPlot$DA <- "Up-regulated"
toPlot$DA[toPlot$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])] <- "Down-regulated"

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_MM_DA.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho, color = factor(DA))) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  scale_color_manual(values = c("#4475B3", "#D7342A")) +
  guides(color = guide_legend(title = "Differential accessibility")) +
  # legend title
  theme_minimal()
dev.off()



toPlot <- data_sig[data_sig$peak_ID%in%SMM_HC,]
toPlot <- as.data.frame(toPlot)
toPlot$DA <- "Up-regulated"
toPlot$DA[toPlot$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_SMM_HC==-1,])] <- "Down-regulated"

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_SMM_DA.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho, color = factor(DA))) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  scale_color_manual(values = c("#4475B3", "#D7342A")) +
  guides(color = guide_legend(title = "Differential accessibility")) +
  # legend title
  theme_minimal()
dev.off()



toPlot <- data_sig[data_sig$peak_ID%in%MGUS_HC,]
toPlot <- as.data.frame(toPlot)
toPlot$DA <- "Up-regulated"
toPlot$DA[toPlot$peak_ID%in%rownames(dea_res_noscore[dea_res_noscore$sig_MGUS_HC==-1,])] <- "Down-regulated"

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_MGUS_DA.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho, color = factor(DA))) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  scale_color_manual(values = c("#4475B3", "#D7342A")) +
  guides(color = guide_legend(title = "Differential accessibility")) +
  # legend title
  theme_minimal()
dev.off()



# Plot non-significant peaks


toPlot <- data_sig[!data_sig$peak_ID%in%unique(c(MGUS_HC,SMM_HC,MM_HC)),]
toPlot <- as.data.frame(toPlot)

library(ggplot2)
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
# plot the corrplot between occupancy and correlation
pdf("occupancy_corrplot_nonDA.pdf", width = 8, height = 5)
ggplot(toPlot, aes(x = occupancy, y = spearman_rho)) +
  geom_point() +
  labs(title = NULL,
       x = "Occupancy (%)",
       y = "Spearman correlation") +
  # legend title
  theme_minimal()
dev.off()


# Save the significant peaks for further analysis
MGUS_HC <- unique(data_sig$peak_ID[data_sig$peak_ID%in%MGUS_HC])
SMM_HC <- unique(data_sig$peak_ID[data_sig$peak_ID%in%SMM_HC])
MM_HC <- unique(data_sig$peak_ID[data_sig$peak_ID%in%MM_HC])

save(MGUS_HC, SMM_HC, MM_HC, file = "coordinated_peaks.RData")




# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
anno <- readRDS("peak_annotation.RDS")
rownames(anno) <- paste0(anno$seqnames, "_", anno$start, "_", anno$end)

# HC -> MGUS
MGUS_HC <- anno[MGUS_HC,]

SMM_HC <- anno[SMM_HC,]

MM_HC <- anno[MM_HC,]



# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data <- rbind(anno, MGUS_HC, SMM_HC, MM_HC)

# add a column with the name of the dataset
data$dataset <- rep(c("ATAC-seq", "MGUS vs HC", "SMM vs HC", "MM vs HC"), c(nrow(anno), nrow(MGUS_HC), nrow(SMM_HC), nrow(MM_HC)))
data$dataset <- factor(data$dataset, levels=c("ATAC-seq", "MGUS vs HC", "SMM vs HC", "MM vs HC"))
# library
library(ggplot2)
library(reshape2)
library(viridis)

library(dplyr)

# mutate data to calculate the percentage of each annotation
data <- data %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
# plot %
pdf("annotation_barplot.pdf", width=5, height=5)
ggplot(data, aes(fill=annotation2, y=freq, x=`dataset`)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T) +
    # legend title
    labs(fill="Region") +
    ylab("Percentage") +
    xlab("") + 
    ggtitle("Peaks' annotation") +
    # change x axis labels 45 degress
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()