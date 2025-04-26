# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")


setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
peaks_association <- readRDS("peaks_association.RDS")

peaks_association$logFC_MGUS_HC <- dea_results[peaks_association$peak_ID,]$logFC_MGUS_HC
peaks_association$logFC_SMM_HC <- dea_results[peaks_association$peak_ID,]$logFC_SMM_HC
peaks_association$logFC_MM_HC <- dea_res_noscore[peaks_association$peak_ID,]$logFC_MM_HC


dim(peaks_association)
length(unique(peaks_association$peak_ID)) # 77990

# Keep only with p.adj < 0.05
peaks_association <- peaks_association[peaks_association$spearman_adj.pval < 0.05,]

dim(peaks_association) # 60358
length(unique(peaks_association$peak_ID)) # 60177

# Plot correlation of spearman rho and logFC for each contrast
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations/Exploration")
library(ggplot2)

pdf("rho_logFC_MGUS.pdf", width = 5, height = 5)
ggplot(peaks_association, aes(x = spearman_rho, y = logFC_MGUS_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MGUS vs HC")
dev.off()


pdf("rho_logFC_SMM.pdf", width = 5, height = 5)
ggplot(peaks_association, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC SMM vs HC")
dev.off()


pdf("rho_logFC_MM.pdf", width = 5, height = 5)
ggplot(peaks_association, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MM vs HC")
dev.off()


# Load ATAC final peaks filtered by global
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")

peaks_association$sig_MGUS_HC <- "ND"
peaks_association[peaks_association$peak_ID%in%MGUS_HC_filter,]$sig_MGUS_HC <- "DA"

peaks_association$sig_SMM_HC <- "ND"
peaks_association[peaks_association$peak_ID%in%SMM_HC_filter,]$sig_SMM_HC <- "DA"

peaks_association$sig_MM_HC <- "ND"
peaks_association[peaks_association$peak_ID%in%MM_HC_filter,]$sig_MM_HC <- "DA"


setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations/Exploration")
toPlot <- peaks_association[peaks_association$sig_MGUS_HC == "DA",]
pdf("rho_logFC_MGUS_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_MGUS_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MGUS vs HC")
dev.off()

toPlot <- peaks_association[peaks_association$sig_SMM_HC == "DA",]
pdf("rho_logFC_SMM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC SMM vs HC")
dev.off()

toPlot <- peaks_association[peaks_association$sig_MM_HC == "DA",]
pdf("rho_logFC_MM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MM vs HC")
dev.off()




# PLOT DA in all contrasts


setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations/Exploration/MGUS_SMM_MM")
toPlot <- peaks_association[peaks_association$sig_MGUS_HC == "DA" & peaks_association$sig_SMM_HC == "DA" & peaks_association$sig_MM_HC == "DA",]
dim(toPlot) # 4621
length(unique(toPlot$peak_ID)) # 4610

pdf("rho_logFC_MGUS_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_MGUS_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MGUS vs HC")
dev.off()

pdf("rho_logFC_SMM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC SMM vs HC")
dev.off()

pdf("rho_logFC_MM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MM vs HC")
dev.off()



setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations/Exploration/SMM_MM")
toPlot <- peaks_association[peaks_association$sig_MGUS_HC == "ND" & peaks_association$sig_SMM_HC == "DA" & peaks_association$sig_MM_HC == "DA",]
dim(toPlot) # 3415
length(unique(toPlot$peak_ID)) # 3404

pdf("rho_logFC_MGUS_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_MGUS_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MGUS vs HC")
dev.off()

pdf("rho_logFC_SMM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC SMM vs HC")
dev.off()

pdf("rho_logFC_MM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MM vs HC")
dev.off()




setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations/Exploration/MM")
toPlot <- peaks_association[peaks_association$sig_MGUS_HC == "ND" & peaks_association$sig_SMM_HC == "ND" & peaks_association$sig_MM_HC == "DA",]
dim(toPlot) # 5479
length(unique(toPlot$peak_ID)) # 5454

pdf("rho_logFC_MGUS_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_MGUS_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MGUS vs HC")
dev.off()

pdf("rho_logFC_SMM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC SMM vs HC")
dev.off()

pdf("rho_logFC_MM_DA.pdf", width = 5, height = 5)
ggplot(toPlot, aes(x = spearman_rho, y = logFC_SMM_HC)) + 
geom_point(alpha = 0.5) + 
theme_minimal() + 
ggtitle("Spearman rho vs logFC MM vs HC")
dev.off()

