####################
### sankey plots ###
####################

#### LOAD THE DATA ####
# Load the DEA results

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
load("final_genes.RData")

genes_toPlot <- unique(c(MGUS_HC_filter, SMM_HC_filter, MM_HC_filter))

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
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/PLOTS")
png("rna_dea_sanky.png", width=1500, height=1000, res=300)
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