# Get peaks from the TF analysis with motifs binding
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
# load pathways
library(fgsea)
library(ggplot2)
library(ggalluvial)

custom_pathways <- gmtPathways("custom_pathways.gmt")

#get TF names
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
TF_data <- read.table("TFs.txt", header = T)
rownames(TF_data) <- make.unique(TF_data$TFs)

# Plot distribution of peaks
# Load the DEA results

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
atac_dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
atac_dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")


# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
rna_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
rna_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load OCR-gene pairs
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_pairs.RData")
OCR_gene_final <- rbind(MGUS_HC, SMM_HC, MM_HC)
# remove duplicated rows
OCR_gene_final <- OCR_gene_final[!duplicated(OCR_gene_final), ]

OCR_gene_final$MGUS_HC_atac <- 0
OCR_gene_final$SMM_HC_atac <- 0
OCR_gene_final$MM_HC_atac <- 0
OCR_gene_final$MGUS_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_dea_results)] <- atac_dea_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_dea_results)], "sig_MGUS_HC"]
OCR_gene_final$SMM_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_dea_results)] <- atac_dea_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_dea_results)], "sig_SMM_HC"]
OCR_gene_final$MM_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_dea_res_noscore)] <- atac_dea_res_noscore[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_dea_res_noscore)], "sig_MM_HC"]

OCR_gene_final$MGUS_HC_rna <- 0
OCR_gene_final$SMM_HC_rna <- 0
OCR_gene_final$MM_HC_rna <- 0
OCR_gene_final$MGUS_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_results)] <- rna_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_results)], "sig_MGUS_HC"]
OCR_gene_final$SMM_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_results)] <- rna_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_results)], "sig_SMM_HC"]
OCR_gene_final$MM_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_res_noscore)] <- rna_res_noscore[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_res_noscore)], "sig_MM_HC"]

# Function sankey
sankey_plot <- function(toPlot, TF) {
    # Reshape for ggalluvial
    genes_long <- tidyr::pivot_longer(
    toPlot,
    cols = c(MGUS_HC_atac, MGUS_HC_rna, SMM_HC_atac, SMM_HC_rna, MM_HC_atac, MM_HC_rna),
    names_to = "Condition",
    values_to = "Category")

    genes_long$Group <- gsub("_(atac|rna)", "", genes_long$Condition)
    genes_long$Modality <- gsub(".*_(atac|rna)", "\\1", genes_long$Condition)

    genes_long$Group <- factor(genes_long$Group, levels = c("MGUS_HC", "SMM_HC", "MM_HC"))
    genes_long$Modality <- factor(genes_long$Modality, levels = c("atac", "rna"))

    genes_long <- genes_long[order(genes_long$Group, genes_long$Modality), ]
    genes_long$Condition <- interaction(genes_long$Group, genes_long$Modality, sep = "_")

    genes_long$Condition <- factor(genes_long$Condition, levels = unique(genes_long$Condition))

    # Add spacer levels
    genes_long$Display_Condition <- as.character(genes_long$Condition)
    genes_long$Display_Condition <- factor(genes_long$Display_Condition,
    levels = c(
    "MGUS_HC_atac", "MGUS_HC_rna", "spacer1",
    "SMM_HC_atac", "SMM_HC_rna", "spacer2",
    "MM_HC_atac", "MM_HC_rna"))

  # Plot
  p <- ggplot2::ggplot(genes_long[!grepl("spacer", genes_long$Display_Condition), ], aes(
  x = Display_Condition, stratum = Category, alluvium = Pair,
  fill = Category, label = Category)) +
  geom_flow(stat = "alluvium", alpha = 0.8, aes(order = as.integer(Display_Condition))) +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "OCR-Gene Pairs Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", 
                               "Down-regulated" = "#4475B3", 
                               "Non-differential" = "#F0F0F0")) +
  scale_x_discrete(drop = FALSE)
  ggsave(p, filename = paste0(TF, "_targets.png"), width = 10, height = 6, dpi = 300)
}

# Function
density_plot_2 <- function(toPlot, contrast) {
    p <- ggplot2::ggplot(toPlot, aes(x = MGUS_HC, fill = "MGUS_HC")) +
    geom_density(alpha = 0.8) +
    geom_density(aes(x = MM_HC, fill = "MM_HC"), alpha = 0.8) +
    geom_density(aes(x = SMM_HC, fill = "SMM_HC"), alpha = 0.8) +
    xlab("log2FC") + 
    ylab("Density") +
    # dash line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("MGUS_HC" = "#D27D46", "SMM_HC" = "#D92F5E", "MM_HC" = "#8D4E85"), labels = c("MGUS vs HC","SMM vs HC", "MM vs HC")) +
    guides(fill = guide_legend(title = "Transition", override.aes = aes(label = ""))) +
    labs(color = "") +
    ggtitle("Density of log2FC") +
    theme_minimal()
    
    ggsave(p, filename = paste0(contrast, "_density_logFC.png"), width = 10, height = 10, dpi = 300)
}

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/05_targets")

TF <- c("JUND.1", "FOS")
for (i in 1:length(TF)) {
  # Get JASPAR ID
  TF_ID <- TF_data[TF[i],]$JASPAR_ID

  # Get peaks with motifs
  TF_motifs <- custom_pathways[[TF_ID]]
  cat("Nr of peaks with", TF[i], "motif:", length(TF_motifs), "\n")

  # Plot distribution of peaks
  toPlot <- data.frame(atac_dea_results[TF_motifs,]$logFC_MGUS_HC,
                        atac_dea_results[TF_motifs,]$logFC_SMM_HC,
                        atac_dea_res_noscore[TF_motifs,]$logFC_MM_HC)
  
  colnames(toPlot) <- c("MGUS_HC", "SMM_HC", "MM_HC")
  rownames(toPlot) <- TF_motifs

  density_plot_2(toPlot, TF[i])

  # Get targets 
  TF_MGUS <- MGUS_HC[MGUS_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_SMM <- SMM_HC[SMM_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_MM <- MM_HC[MM_HC$peak_ID %in% TF_motifs,]$gene_ID 

  cat(TF[i], "targets HCvsMGUS:", length(TF_MGUS), "\n")
  cat(TF[i], "targets HCvsSMM:", length(TF_SMM), "\n")
  cat(TF[i], "targets HCvsMM:", length(TF_MM), "\n")

  # Plot sankey 
  # OCR-gene pairs with peaks that have motif for the TF + targets 
  #### PREP THE DATA ####
  dynamics <- OCR_gene_final[OCR_gene_final$peak_ID%in%TF_motifs,]
  rownames(dynamics) <- dynamics$pairs
  toPlot <- dynamics[,c("MGUS_HC_atac", "MGUS_HC_rna", "SMM_HC_atac", "SMM_HC_rna", "MM_HC_atac", "MM_HC_rna")]

  # Change 0s to "non differential"
  toPlot[toPlot == 0] <- "Non-differential"
  #Change -1 to "Down"
  toPlot[toPlot == -1] <- "Down-regulated"
  #Change 1 to "Up"
  toPlot[toPlot == 1] <- "Up-regulated"

  toPlot$Pair <- rownames(dynamics)

  sankey_plot(toPlot, TF[i])
}


# Check JUND::FOS motifs in the same peaks


TF_jund <- TF_data["JUND.1",]$JASPAR_ID
TF_fos <- TF_data["FOS",]$JASPAR_ID

# Get peaks with motifs
TF_motifs_jund <- custom_pathways[[TF_jund]]
TF_motifs_fos <- custom_pathways[[TF_fos]]

common <- intersect(TF_motifs_jund, TF_motifs_fos) #192

# Get targets 
TF_MGUS <- MGUS_HC[MGUS_HC$peak_ID %in% common,]$gene_ID 
TF_SMM <- SMM_HC[SMM_HC$peak_ID %in% common,]$gene_ID 
TF_MM <- MM_HC[MM_HC$peak_ID %in% common,]$gene_ID 

cat("JUND::FOS targets HCvsMGUS:", length(TF_MGUS), "\n") #26
cat("JUND::FOS HCvsSMM:", length(TF_SMM), "\n") #61
cat("JUND::FOS HCvsMM:", length(TF_MM), "\n") #221

together <- unique(c(TF_MGUS, TF_SMM, TF_MM)) # 208
# Plot sankey 
# OCR-gene pairs with peaks that have motif for the TF + targets 
#### PREP THE DATA ####
dynamics <- OCR_gene_final[OCR_gene_final$peak_ID%in%common,]
rownames(dynamics) <- dynamics$pairs
toPlot <- dynamics[,c("MGUS_HC_atac", "MGUS_HC_rna", "SMM_HC_atac", "SMM_HC_rna", "MM_HC_atac", "MM_HC_rna")]

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Pair <- rownames(dynamics)

sankey_plot(toPlot, "JUND::FOS")

table(TF_MGUS%in%TF_SMM) # 22 yes
table(TF_MGUS%in%TF_MM) # 23 yes
table(TF_SMM%in%TF_MM) # 52 yes





  TF_ID <- TF_data["JUND.1",]$JASPAR_ID

  # Get peaks with motifs
  TF_motifs <- custom_pathways[[TF_ID]]

  # Get targets 
  TF_MGUS <- MGUS_HC[MGUS_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_SMM <- SMM_HC[SMM_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_MM <- MM_HC[MM_HC$peak_ID %in% TF_motifs,]$gene_ID 

jund_targets <- unique(c(TF_MGUS, TF_SMM, TF_MM))

  TF_ID <- TF_data["FOS",]$JASPAR_ID

  # Get peaks with motifs
  TF_motifs <- custom_pathways[[TF_ID]]

  # Get targets 
  TF_MGUS <- MGUS_HC[MGUS_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_SMM <- SMM_HC[SMM_HC$peak_ID %in% TF_motifs,]$gene_ID 
  TF_MM <- MM_HC[MM_HC$peak_ID %in% TF_motifs,]$gene_ID 

fos_targets <- unique(c(TF_MGUS, TF_SMM, TF_MM))

table(intersect(jund_targets, fos_targets)%in%together)