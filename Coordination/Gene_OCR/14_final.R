### Check if annotated genes in differential peaks are differentially expressed 

rm(list = ls())
#################
# LOAD RNA DATA #
#################
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
load("final_genes.RData")
MGUS_HC_genes <- MGUS_HC_filter
SMM_HC_genes <- SMM_HC_filter
MM_HC_genes <- MM_HC_filter

##################
# LOAD ATAC DATA #
##################
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
peak_annotation <- readRDS("peak_annotation.RDS")

# DEA RESULTS - filtered by ChIP
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")
MGUS_HC_atac <- MGUS_HC
SMM_HC_atac <- SMM_HC
MM_HC_atac <- MM_HC

#################
### PROMOTERS ###
#################

# For promoter OCRs, we directly annotated OCRs to genes based on TSS distance (<1kb)

promoters <- peak_annotation[peak_annotation$annotation2 == "Promoter",]
dim(promoters) # 14747 OCRs mapping in promoter regions
# remove rows with NA genes
promoters <- promoters[!is.na(promoters$ENSEMBL),]
# How many OCR-gene pairs ?
dim(promoters) # 13127 OCRs mapping in promoter regions 
summary(is.na(promoters$name))

#################
### ENHANCERS ###
#################

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/GRN/General/enhancers")
enhancers <- read.table("OCR_gene_association.txt", header = TRUE, sep = "\t")
# keep only significant associations p.adj < 0.05
enhancers <- enhancers[enhancers$spearman_adj.pval < 0.05,]
dim(enhancers) # # 39421 OCRs mapping in enhancer regions 
enhancers <- enhancers[!is.na(enhancers$gene_ID),]
# How many OCR-gene pairs ?
dim(enhancers) # # 38549 OCRs mapping in enhancer regions 
summary(is.na(enhancers$peak_ID))

# Are there any peaks in promoters and enhancers
summary(enhancers$peak_ID %in% promoters$name) # NO


# Create OCR_gene pairs for promoters and enhancers
promoters <- promoters[,c("name", "ENSEMBL")]
colnames(promoters) <- c("peak_ID", "gene_ID")
enhancers <- enhancers[,c("peak_ID", "gene_ID")]

# Merge OCR_gene pairs
OCR_gene <- rbind(promoters, enhancers)
dim(OCR_gene) # 51676 OCR-gene pairs

OCR_gene$pairs <- paste0(OCR_gene$peak_ID, "-", OCR_gene$gene_ID)
summary(duplicated(OCR_gene$pairs)) # NO

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
saveRDS(OCR_gene, file = "OCR_gene_pairs.RDS")
# OCR_gene pairs for differential peaks
DA_MGUS_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% MGUS_HC_atac,]
dim(DA_MGUS_OCR_gene) # 1922 OCR-gene pairs

DA_SMM_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% SMM_HC_atac,]
dim(DA_SMM_OCR_gene) # 3183 OCR-gene pairs

DA_MM_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% MM_HC_atac,]
dim(DA_MM_OCR_gene) # 5676 OCR-gene pairs

# save 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
saveRDS(DA_MGUS_OCR_gene, file = "DA_MGUS_OCR_gene.RDS")
saveRDS(DA_SMM_OCR_gene, file = "DA_SMM_OCR_gene.RDS")
saveRDS(DA_MM_OCR_gene, file = "DA_MM_OCR_gene.RDS")

# OCR_gene pairs for differential peaks and differential genes
DA_MGUS_OCR_gene_DE <- DA_MGUS_OCR_gene[DA_MGUS_OCR_gene$gene_ID %in% MGUS_HC_genes,]
dim(DA_MGUS_OCR_gene_DE) # 246 OCR-gene pairs

DA_SMM_OCR_gene_DE <- DA_SMM_OCR_gene[DA_SMM_OCR_gene$gene_ID %in% SMM_HC_genes,]
dim(DA_SMM_OCR_gene_DE) # 617 OCR-gene pairs

DA_MM_OCR_gene_DE <- DA_MM_OCR_gene[DA_MM_OCR_gene$gene_ID %in% MM_HC_genes,]
dim(DA_MM_OCR_gene_DE) # 1879 OCR-gene pairs


library(ggVennDiagram)
library(ggplot2)
### FILTER DE genes if it is in a pair with DA peak
# Create a list of gene sets
genes_stes <- list(MGUS_HC_genes, unique(DA_MGUS_OCR_gene_DE$gene_ID))

# Create a list of gene set names
genes_names <- c("DE genes","Genes from OCR-gene pairs")


# Save final intersections
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
MGUS_HC <- intersect(MGUS_HC_genes, unique(DA_MGUS_OCR_gene_DE$gene_ID))
SMM_HC <- intersect(SMM_HC_genes, unique(DA_SMM_OCR_gene_DE$gene_ID))
MM_HC <- intersect(MM_HC_genes, unique(DA_MM_OCR_gene_DE$gene_ID))

save(MGUS_HC, SMM_HC, MM_HC, file = "coordination_genes.RData")


library(ggVennDiagram)
library(ggplot2)
# Create a list of gene sets
genes_stes <- list(MGUS_HC, SMM_HC, MM_HC)

# Create a list of gene set names
genes_names <- c("MGUS vs HC", "SMM vs HC", "MM vs HC")

# Create a Venn diagram
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
png("venn_diagram_RNA.png", width=2000, height=2000, res=300)
ggVennDiagram(genes_stes, genes_names) + 
# fill gradient based on % of overlap values min = 0, max = 100%
scale_fill_gradientn(colours = c("white", "#CC6633", "#CC3333"), values = scales::rescale(c(0, mean(length(MM_HC), length(MM_HC)))) ) 
dev.off()



# Save final intersections
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
MGUS_HC <- intersect(MGUS_HC_atac, unique(DA_MGUS_OCR_gene_DE$peak_ID))
SMM_HC <- intersect(SMM_HC_atac, unique(DA_SMM_OCR_gene_DE$peak_ID))
MM_HC <- intersect(MM_HC_atac, unique(DA_MM_OCR_gene_DE$peak_ID))

save(MGUS_HC, SMM_HC, MM_HC, file = "coordination_peaks.RData")


library(ggVennDiagram)
library(ggplot2)
# Create a list of gene sets
genes_stes <- list(MGUS_HC, SMM_HC, MM_HC)

# Create a list of gene set names
genes_names <- c("MGUS vs HC", "SMM vs HC", "MM vs HC")

# Create a Venn diagram
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
png("venn_diagram_ATAC.png", width=2000, height=2000, res=300)
ggVennDiagram(genes_stes, genes_names) + 
# fill gradient based on % of overlap values min = 0, max = 100%
scale_fill_gradientn(colours = c("white", "#CC6633", "#CC3333"), values = scales::rescale(c(0, mean(length(MM_HC), length(MM_HC)))) ) 
dev.off()



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

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
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





# Save final pairs
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
MGUS_HC <- DA_MGUS_OCR_gene_DE
SMM_HC <- DA_SMM_OCR_gene_DE
MM_HC <- DA_MM_OCR_gene_DE

# save final pairs
save(MGUS_HC, SMM_HC, MM_HC, file = "coordination_pairs.RData")


# Load the DEA results
##### ATAC
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_atac <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_noscore_atac <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

##### RNA
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_rna <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_noscore_rna <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

MM_HC$atac_logFC <- dea_noscore_atac[MM_HC[,1],]$logFC_MM_HC
MM_HC$rna_logFC <- dea_noscore_rna[MM_HC[,2],]$logFC_MM_HC

SMM_HC$atac_logFC <- dea_atac[SMM_HC[,1],]$logFC_SMM_HC
SMM_HC$rna_logFC <- dea_rna[SMM_HC[,2],]$logFC_SMM_HC

MGUS_HC$atac_logFC <- dea_atac[MGUS_HC[,1],]$logFC_MGUS_HC
MGUS_HC$rna_logFC <- dea_rna[MGUS_HC[,2],]$logFC_MGUS_HC

# Correlation plots of logFCs
library(ggplot2)

# Function
cor_plot <- function(toPlot, contrast, annotation) {
    # keep only promoters or enhancers 
    if ( annotation == "promoters") {
        data <- toPlot[toPlot$peak_ID %in% promoters$peak_ID,]
    } else {
        data <- toPlot[toPlot$peak_ID %in% enhancers$peak_ID,]
    }
    data <- data.frame(data$atac_logFC, data$rna_logFC)
    rownames(data) <- 1:nrow(data)
    colnames(data) <- c("OCRs", "Genes")

    correlation_coefficient <- cor(data$OCRs, data$Genes, method = "spearman")

    p <-  ggplot(data, aes(x = Genes, y = OCRs)) +
            geom_point() +  # Plot points
            # add line in x = 0 
            geom_vline(xintercept = 0, color = "black") +
            geom_hline(yintercept = 0, color = "black") +
            annotate("text", x = min(data$Genes)+0.5, y = max(data$OCRs), label = paste("rho=", round(correlation_coefficient, 5))) +
            xlab("Genes") + 
            ylab(paste0("OCRs - ", annotation)) +
            ggtitle(paste0("Correlation of log2FC - ", contrast)) +
            theme_minimal()
    ggsave(p, filename = paste0("correlation_",contrast, "_", annotation, ".png"), width = 10, height = 10, dpi = 300)
    }





# correlation plot
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/Correlations")
cor_plot(MGUS_HC, "MGUS_HC", "enhancers")
cor_plot(MGUS_HC, "MGUS_HC", "promoters")
cor_plot(SMM_HC, "SMM_HC", "enhancers")
cor_plot(SMM_HC, "SMM_HC", "promoters")
cor_plot(MM_HC, "MM_HC", "enhancers")
cor_plot(MM_HC, "MM_HC", "promoters")

