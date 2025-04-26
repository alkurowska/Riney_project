### Check if annotated genes in differential peaks are differentially expressed 

rm(list = ls())
#################
# LOAD RNA DATA #
#################
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
mouse_dea <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
dim(mouse_dea) # 19423

dynamics <- mouse_dea[,c("sig_CyclinD1_MGUS", "sig_Mmset_MGUS", "sig_Trp53_MGUS",  "sig_MIc_MGUS", 
                        "sig_CyclinD1_MM", "sig_Mmset_MM", "sig_Trp53_MM", "sig_MIc_MM")]

# Remove all zeros
dynamics <- dynamics[rowSums(dynamics) != 0,]
dim(dynamics) # 9229

# PLOT sankey

#### PREP THE DATA ####
toPlot <- dynamics

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Gene <- rownames(dynamics)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial

genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(sig_CyclinD1_MGUS, sig_Mmset_MGUS, sig_Trp53_MGUS, sig_MIc_MGUS,
           sig_CyclinD1_MM, sig_Mmset_MM, sig_Trp53_MM, sig_MIc_MM),
  names_to = "Condition",
  values_to = "Category"
)

genes_long$Condition <- gsub("sig_", "", genes_long$Condition)

genes_long$Group <- "slow"
genes_long$Group[grepl("Trp53|MIc", genes_long$Condition)] <- "rapid"

genes_long$Condition <- gsub("MGUS", "MGUS_HC", genes_long$Condition)
genes_long$Condition <- gsub("MM", "MM_HC", genes_long$Condition)

genes_long$Contrast <- gsub(".*_(MGUS_HC|MM_HC)", "\\1", genes_long$Condition)
genes_long$Contrast <- factor(genes_long$Contrast, levels = c("MGUS_HC","MM_HC"))

# Add spacer levels
genes_long$Display_Condition <- as.character(genes_long$Condition)
genes_long$Display_Condition <- factor(genes_long$Display_Condition,
  levels = c(
    "CyclinD1_MGUS_HC", "Mmset_MGUS_HC", 
    "Trp53_MGUS_HC",  "MIc_MGUS_HC", "spacer1",
    "CyclinD1_MM_HC", "Mmset_MM_HC", 
    "Trp53_MM_HC", "MIc_MM_HC"
  )
)

# Plot

png("gene_sanky.png", width=2000, height=1000, res=300)

# Plot with spacer strata dropped
ggplot(genes_long[!grepl("spacer", genes_long$Display_Condition), ], aes(
  x = Display_Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium", alpha = 0.8, aes(order = as.integer(Display_Condition))) +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", 
                               "Down-regulated" = "#4475B3", 
                               "Non-differential" = "#F0F0F0")) +
  scale_x_discrete(drop = FALSE)
dev.off()



# Plot per model
models <- c("CyclinD1", "Mmset", "Trp53", "MIc")

for(i in 1:length(models)){
  model <- models[i]
  
  # Filter for model
  toPlot <- dynamics[,which(grepl(model, colnames(dynamics)))]
  
  # Change 0s to "non differential"
  toPlot[toPlot == 0] <- "Non-differential"
  #Change -1 to "Down"
  toPlot[toPlot == -1] <- "Down-regulated"
  #Change 1 to "Up"
  toPlot[toPlot == 1] <- "Up-regulated"
  
  toPlot$Gene <- rownames(dynamics)
  
  # Reshape for ggalluvial

  genes_long <- tidyr::pivot_longer(
   toPlot,
   cols = colnames(toPlot)[1:2],
   names_to = "Condition",
   values_to = "Category"
  )

  genes_long$Condition <- gsub("sig_", "", genes_long$Condition)

  genes_long$Condition <- gsub("MGUS", "MGUS_HC", genes_long$Condition)
  genes_long$Condition <- gsub("MM", "MM_HC", genes_long$Condition)

  # Plot
  p <- ggplot(genes_long, aes(
  x = Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium", alpha = 0.8) +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", 
                               "Down-regulated" = "#4475B3", 
                               "Non-differential" = "#F0F0F0")) +
  scale_x_discrete(drop = FALSE)
ggsave(paste0("gene_sanky_", model, ".png"), plot = p, width=7, height=5, dpi=300)
}

