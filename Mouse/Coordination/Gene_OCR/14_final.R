### Check if annotated genes in differential peaks are differentially expressed 

rm(list = ls())
#################
# LOAD RNA DATA #
#################
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

CyclinD1_MGUS_g  <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS_g  <- rownames(dea_results[dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS_g  <- rownames(dea_results[dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS_g <- rownames(dea_results[dea_results$sig_MIc_MGUS != 0,])

genes_MGUS <- unique(c(CyclinD1_MGUS_g, Mmset_MGUS_g, Trp53_MGUS_g, MIc_MGUS_g))

CyclinD1_MM_g  <- rownames(dea_results[dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM_g  <- rownames(dea_results[dea_results$sig_Mmset_MM != 0,])
Trp53_MM_g  <- rownames(dea_results[dea_results$sig_Trp53_MM != 0,])
MIc_MM_g <- rownames(dea_results[dea_results$sig_MIc_MM != 0,])

genes_MM <- unique(c(CyclinD1_MM_g, Mmset_MM_g, Trp53_MM_g, MIc_MM_g))

##################
# LOAD ATAC DATA #
##################
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Peaks")
peak_annotation <- readRDS("peak_annotation.RDS")

# DEA RESULTS - filtered by ChIP
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA")
atac_dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

CyclinD1_MGUS  <- rownames(atac_dea_results[atac_dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS  <- rownames(atac_dea_results[atac_dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS  <- rownames(atac_dea_results[atac_dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS <- rownames(atac_dea_results[atac_dea_results$sig_MIc_MGUS != 0,])

all_MGUS <- unique(c(CyclinD1_MGUS, Mmset_MGUS, Trp53_MGUS, MIc_MGUS))

CyclinD1_MM  <- rownames(atac_dea_results[atac_dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM  <- rownames(atac_dea_results[atac_dea_results$sig_Mmset_MM != 0,])
Trp53_MM  <- rownames(atac_dea_results[atac_dea_results$sig_Trp53_MM != 0,])
MIc_MM <- rownames(atac_dea_results[atac_dea_results$sig_MIc_MM != 0,])

all_MM <- unique(c(CyclinD1_MM, Mmset_MM, Trp53_MM, MIc_MM))

#################
### PROMOTERS ###
#################

# For promoter OCRs, we directly annotated OCRs to genes based on TSS distance (<1kb)

promoters <- peak_annotation[peak_annotation$annotation2 == "Promoter",]
dim(promoters) # 11963 OCRs mapping in promoter regions
# remove rows with NA genes
promoters <- promoters[!is.na(promoters$ENSEMBL),]
# How many OCR-gene pairs ?
dim(promoters) # 11646 OCRs mapping in promoter regions 
summary(is.na(promoters$name))

#################
### ENHANCERS ###
#################

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/GRN/General/enhancers")
enhancers <- read.table("OCR_gene_association.txt", header = TRUE, sep = "\t")
# keep only significant associations p.adj < 0.05
enhancers <- enhancers[enhancers$spearman_adj.pval < 0.05,]
dim(enhancers) # # 14777 OCRs mapping in enhancer regions 
enhancers <- enhancers[!is.na(enhancers$gene_ID),]
# How many OCR-gene pairs ?
dim(enhancers) # # 13878 OCRs mapping in enhancer regions 
summary(is.na(enhancers$peak_ID))

# Are there any peaks in promoters and enhancers
summary(enhancers$peak_ID %in% promoters$name) # NO


# Create OCR_gene pairs for promoters and enhancers
promoters <- promoters[,c("name", "ENSEMBL")]
colnames(promoters) <- c("peak_ID", "gene_ID")
enhancers <- enhancers[,c("peak_ID", "gene_ID")]

# Merge OCR_gene pairs
OCR_gene <- rbind(promoters, enhancers)
dim(OCR_gene) # 25524 OCR-gene pairs

OCR_gene$pairs <- paste0(OCR_gene$peak_ID, "-", OCR_gene$gene_ID)
summary(duplicated(OCR_gene$pairs)) # NO

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
saveRDS(OCR_gene, file = "OCR_gene_pairs.RDS")
# OCR_gene pairs for differential peaks
DA_MGUS_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% all_MGUS,]
dim(DA_MGUS_OCR_gene) # 1799 OCR-gene pairs

DA_MM_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% all_MM,]
dim(DA_MM_OCR_gene) # 2522 OCR-gene pairs


# OCR_gene pairs for differential peaks and differential genes
DA_MGUS_OCR_gene_DE <- DA_MGUS_OCR_gene[DA_MGUS_OCR_gene$gene_ID %in% genes_MGUS,]
dim(DA_MGUS_OCR_gene_DE) # 590 OCR-gene pairs

DA_MM_OCR_gene_DE <- DA_MM_OCR_gene[DA_MM_OCR_gene$gene_ID %in% genes_MM,]
dim(DA_MM_OCR_gene_DE) # 1676 OCR-gene pairs


# Save final intersections
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
MGUS_HC <- intersect(all_MGUS, unique(DA_MGUS_OCR_gene_DE$peak_ID))
MM_HC <- intersect(all_MM, unique(DA_MM_OCR_gene_DE$peak_ID))

save(MGUS_HC, MM_HC, file = "coordination_peaks.RData")

MGUS_HC <- intersect(genes_MGUS, unique(DA_MGUS_OCR_gene_DE$gene_ID))
MM_HC <- intersect(genes_MM, unique(DA_MM_OCR_gene_DE$gene_ID))

save(MGUS_HC, MM_HC, file = "coordination_genes.RData")

# Save final pairs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
MGUS_HC <- DA_MGUS_OCR_gene_DE
MM_HC <- DA_MM_OCR_gene_DE

# save final pairs
save(MGUS_HC, MM_HC, file = "coordination_pairs.RData")




# MM_HC$atac_logFC <- dea_noscore_atac[MM_HC[,1],]$logFC_MM_HC
# MM_HC$rna_logFC <- dea_noscore_rna[MM_HC[,2],]$logFC_MM_HC

# SMM_HC$atac_logFC <- dea_atac[SMM_HC[,1],]$logFC_SMM_HC
# SMM_HC$rna_logFC <- dea_rna[SMM_HC[,2],]$logFC_SMM_HC

# MGUS_HC$atac_logFC <- dea_atac[MGUS_HC[,1],]$logFC_MGUS_HC
# MGUS_HC$rna_logFC <- dea_rna[MGUS_HC[,2],]$logFC_MGUS_HC

# # Correlation plots of logFCs
# library(ggplot2)

# # Function
# cor_plot <- function(toPlot, contrast, annotation) {
#     # keep only promoters or enhancers 
#     if ( annotation == "promoters") {
#         data <- toPlot[toPlot$peak_ID %in% promoters$peak_ID,]
#     } else {
#         data <- toPlot[toPlot$peak_ID %in% enhancers$peak_ID,]
#     }
#     data <- data.frame(data$atac_logFC, data$rna_logFC)
#     rownames(data) <- 1:nrow(data)
#     colnames(data) <- c("OCRs", "Genes")

#     correlation_coefficient <- cor(data$OCRs, data$Genes, method = "spearman")

#     p <-  ggplot(data, aes(x = Genes, y = OCRs)) +
#             geom_point() +  # Plot points
#             # add line in x = 0 
#             geom_vline(xintercept = 0, color = "black") +
#             geom_hline(yintercept = 0, color = "black") +
#             annotate("text", x = min(data$Genes)+0.5, y = max(data$OCRs), label = paste("rho=", round(correlation_coefficient, 5))) +
#             xlab("Genes") + 
#             ylab(paste0("OCRs - ", annotation)) +
#             ggtitle(paste0("Correlation of log2FC - ", contrast)) +
#             theme_minimal()
#     ggsave(p, filename = paste0("correlation_",contrast, "_", annotation, ".png"), width = 10, height = 10, dpi = 300)
#     }





# # correlation plot
# setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/Correlations")
# cor_plot(MGUS_HC, "MGUS_HC", "enhancers")
# cor_plot(MGUS_HC, "MGUS_HC", "promoters")
# cor_plot(SMM_HC, "SMM_HC", "enhancers")
# cor_plot(SMM_HC, "SMM_HC", "promoters")
# cor_plot(MM_HC, "MM_HC", "enhancers")
# cor_plot(MM_HC, "MM_HC", "promoters")

