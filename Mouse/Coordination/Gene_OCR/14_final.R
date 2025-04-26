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

promoters$name <- paste0(promoters$seqnames, "_", promoters$start, "_", promoters$end)
# Create OCR_gene pairs for promoters and enhancers
promoters <- promoters[,c("name", "ENSEMBL")]
colnames(promoters) <- c("peak_ID", "gene_ID")
enhancers <- enhancers[,c("peak_ID", "gene_ID")]

# Are there any peaks in promoters and enhancers
summary(enhancers$peak_ID %in% promoters$name) # NO


# Merge OCR_gene pairs
OCR_gene <- rbind(promoters, enhancers)
dim(OCR_gene) # 25524 OCR-gene pairs

OCR_gene$pairs <- paste0(OCR_gene$peak_ID, "-", OCR_gene$gene_ID)
summary(duplicated(OCR_gene$pairs)) # NO

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
saveRDS(OCR_gene, file = "OCR_gene_pairs.RDS")
# OCR_gene pairs for differential peaks
DA_MGUS_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% all_MGUS,]
dim(DA_MGUS_OCR_gene) # 870 OCR-gene pairs

DA_MM_OCR_gene <- OCR_gene[OCR_gene$peak_ID %in% all_MM,]
dim(DA_MM_OCR_gene) # 1598 OCR-gene pairs

toSave <- rbind(DA_MGUS_OCR_gene, DA_MM_OCR_gene)
# Remove duplicates
toSave <- toSave[!duplicated(toSave$pairs),]
# SAVE OCR_gene pairs for differential peaks
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
saveRDS(toSave, file = "DA_MGUS_OCR_gene.RDS")

# OCR_gene pairs for differential peaks and differential genes
DA_MGUS_OCR_gene_DE <- DA_MGUS_OCR_gene[DA_MGUS_OCR_gene$gene_ID %in% genes_MGUS,]
dim(DA_MGUS_OCR_gene_DE) # 140 OCR-gene pairs

DA_MM_OCR_gene_DE <- DA_MM_OCR_gene[DA_MM_OCR_gene$gene_ID %in% genes_MM,]
dim(DA_MM_OCR_gene_DE) # 768 OCR-gene pairs


# Save final intersections
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
MGUS_HC <- intersect(all_MGUS, unique(DA_MGUS_OCR_gene_DE$peak_ID))
MM_HC <- intersect(all_MM, unique(DA_MM_OCR_gene_DE$peak_ID))

MGUS_HC <- intersect(genes_MGUS, unique(DA_MGUS_OCR_gene_DE$gene_ID))
MM_HC <- intersect(genes_MM, unique(DA_MM_OCR_gene_DE$gene_ID))

save(MGUS_HC, MM_HC, file = "coordination_genes.RData")

# Save final pairs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
MGUS_HC <- DA_MGUS_OCR_gene_DE
MM_HC <- DA_MM_OCR_gene_DE

# save final pairs
save(MGUS_HC, MM_HC, file = "coordination_pairs.RData")




# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA/")
atac_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
rna_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

OCR_gene_final <- rbind(MGUS_HC, MM_HC)
OCR_gene_final <- OCR_gene_final[!duplicated(OCR_gene_final$pairs),]
dim(OCR_gene_final) # 788

OCR_gene_final$atac_CyclindD1_MGUS <- 0
OCR_gene_final$atac_Mmset_MGUS <- 0
OCR_gene_final$atac_Trp53_MGUS <- 0
OCR_gene_final$atac_MIc_MGUS <- 0
OCR_gene_final$atac_CyclindD1_MM <- 0
OCR_gene_final$atac_Mmset_MM <- 0
OCR_gene_final$atac_Trp53_MM <- 0
OCR_gene_final$atac_MIc_MM <- 0
OCR_gene_final$atac_CyclindD1_MGUS[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_CyclinD1_MGUS"]
OCR_gene_final$atac_Mmset_MGUS[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_Mmset_MGUS"]
OCR_gene_final$atac_Trp53_MGUS[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_Trp53_MGUS"]
OCR_gene_final$atac_MIc_MGUS[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_MIc_MGUS"]
OCR_gene_final$atac_CyclindD1_MM[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_CyclinD1_MM"]
OCR_gene_final$atac_Mmset_MM[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_Mmset_MM"]
OCR_gene_final$atac_Trp53_MM[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_Trp53_MM"]
OCR_gene_final$atac_MIc_MM[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_MIc_MM"]

OCR_gene_final$rna_CyclindD1_MGUS <- 0
OCR_gene_final$rna_Mmset_MGUS <- 0
OCR_gene_final$rna_Trp53_MGUS <- 0
OCR_gene_final$rna_MIc_MGUS <- 0
OCR_gene_final$rna_CyclindD1_MM <- 0
OCR_gene_final$rna_Mmset_MM <- 0
OCR_gene_final$rna_Trp53_MM <- 0
OCR_gene_final$rna_MIc_MM <- 0
OCR_gene_final$rna_CyclindD1_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_CyclinD1_MGUS"]
OCR_gene_final$rna_Mmset_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Mmset_MGUS"]
OCR_gene_final$rna_Trp53_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Trp53_MGUS"]
OCR_gene_final$rna_MIc_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_MIc_MGUS"]
OCR_gene_final$rna_CyclindD1_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_CyclinD1_MM"]
OCR_gene_final$rna_Mmset_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Mmset_MM"]
OCR_gene_final$rna_Trp53_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Trp53_MM"]
OCR_gene_final$rna_MIc_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_MIc_MM"]

head(OCR_gene_final)

library(ggplot2)
library(ggalluvial)

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
# PER MODEL
models <- c("CyclindD1", "Mmset", "Trp53", "MIc")
# summarize the number of genes in each category
for(i in 1:length(models)){
    model <- models[i]
    dynamics <- OCR_gene_final[,grep(model, colnames(OCR_gene_final))]
    dynamics <- dynamics[rowSums(dynamics) != 0,]

    toPlot <- dynamics
    # Change 0s to "non differential"
    toPlot[toPlot == 0] <- "Non-differential"
    #Change -1 to "Down"
    toPlot[toPlot == -1] <- "Down-regulated"
    #Change 1 to "Up"
    toPlot[toPlot == 1] <- "Up-regulated"

    toPlot$Pair <- rownames(dynamics)
    colnames(toPlot) <- gsub(paste0(model, "_"), "", colnames(toPlot))

    # Reshape for ggalluvial

    genes_long <- tidyr::pivot_longer(
     toPlot,
     cols = c(atac_MGUS, atac_MM, rna_MGUS, rna_MM),
     names_to = "Condition",
     values_to = "Category"
    )

    genes_long$Group <- gsub("(atac|rna)_", "", genes_long$Condition)
    genes_long$Modality <- gsub("(atac|rna)_.*", "\\1", genes_long$Condition)

    genes_long$Group <- factor(genes_long$Group, levels = c("MGUS", "MM"))
    genes_long$Modality <- factor(genes_long$Modality, levels = c("atac", "rna"))

    genes_long <- genes_long[order(genes_long$Group, genes_long$Modality), ]
    genes_long$Condition <- interaction(genes_long$Group, genes_long$Modality, sep = "_")

    genes_long$Condition <- factor(genes_long$Condition, levels = unique(genes_long$Condition))

    # Add spacer levels
    genes_long$Display_Condition <- as.character(genes_long$Condition)
    genes_long$Display_Condition <- factor(genes_long$Display_Condition,
    levels = c(
    "MGUS_atac", "MGUS_rna", "spacer1",
    "MM_atac", "MM_rna"))

# Plot

        # Plot with spacer strata dropped
    p <- ggplot(genes_long[!grepl("spacer", genes_long$Display_Condition), ], aes(
     x = Display_Condition, stratum = Category, alluvium = Pair,
     fill = Category, label = Category
        )) +
     geom_flow(stat = "alluvium", alpha = 0.8, aes(order = as.integer(Display_Condition))) +
     geom_stratum() +
    theme_classic() +
    labs(x = NULL, y = "OCR-Gene Pairs Count") +
    scale_fill_manual(values = c("Up-regulated" = "#D7342A", 
                               "Down-regulated" = "#4475B3", 
                               "Non-differential" = "#F0F0F0")) +
    scale_x_discrete(drop = FALSE)
    ggsave(paste0("Sankey_", model, ".pdf"), plot = p, width = 6, height = 4)

}