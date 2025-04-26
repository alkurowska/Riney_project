##### READ DATA #####
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/all_peaks/04_ranking")
MGUS <- read.table("Rank_MGUS.txt",header = T, sep = "\t")
SMM <- read.table("Rank_SMM.txt",header = T, sep = "\t")
MM <- read.table("Rank_MM.txt",header = T, sep = "\t")


# FDR 
MGUS$p.adjust <- p.adjust(MGUS$Score, method = "fdr")
SMM$p.adjust <- p.adjust(SMM$Score, method = "fdr")
MM$p.adjust <- p.adjust(MM$Score, method = "fdr")

# sort by name
MGUS <- MGUS[order(MGUS[,1]),]
SMM <- SMM[order(SMM[,1]),]
MM <- MM[order(MM[,1]),]
summary(rownames(MGUS) == rownames(SMM)) # TRUE
summary(rownames(MGUS) == rownames(MM)) # TRUE

toPlot <- data.frame(
  TF = rownames(MGUS),
  MGUS = MGUS$Score,
  SMM = SMM$Score,
  MM = MM$Score
)

# keep only significant score < 0.05
toPlot <- toPlot[toPlot$MGUS < 0.05 | toPlot$SMM < 0.05 | toPlot$MM < 0.05,]

# Load Pioneers and Crisprs
# List of essential TFs
Pioneer <- as.character( read.table("Pioneers.txt",header = F)[,1])
Crisp <- as.character(read.table("Crispr.txt",header = F)[,1])

# Are any of the essential TFs in the list?
Pioneer[(Pioneer %in% toPlot$TF)]
Crisp[(Crisp %in% toPlot$TF)]


# Plot a box plots
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

color_trans <- c("#D27D46", "#D92F5E", "#8D4E85") 
names(color_trans) <- c("MGUS", "SMM", "MM")

metadata <- data.frame(colnames(toPlot)[2:4])
colnames(metadata) <- c("Transition")

ha1 <- rowAnnotation(
  Transition = metadata$Transition,
  col = list(Transition = as.factor(color_trans)),
  simple_anno_size = unit(1, "cm"),
  na_col = "#D9D9D9", annotation_name_gp= gpar(fontsize = 15),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15, fontface = "bold")),
  show_annotation_name = TRUE)

# Create significance matrix
# Create significance matrix - note the transpose
sig_matrix <- t(matrix(0, nrow = 3, ncol = nrow(toPlot)))
colnames(sig_matrix) <- c("MGUS", "SMM", "MM")
rownames(sig_matrix) <- toPlot$TF
sig_matrix[,1] <- as.numeric(MGUS[toPlot$TF,]$p.adjust < 0.05)
sig_matrix[,2] <- as.numeric(SMM[toPlot$TF,]$p.adjust < 0.05)
sig_matrix[,3] <- as.numeric(MM[toPlot$TF,]$p.adjust < 0.05)

# Function to add dots for significant values
cell_fun = function(j, i, x, y, width, height, fill) {
    if(sig_matrix[j, i] == 1) {
        grid.circle(x = x, y = y, r = unit(1, "mm"), 
                   gp = gpar(col = "black", fill = "black"))
    }
}

# Create dot legend
dot_legend = Legend(
    labels = "FDR < 0.05",
    type = "points",
    pch = 19,
    legend_gp = gpar(col = "black"),
    title = "Significance",
    title_gp = gpar(fontsize = 15, fontface = "bold"),
    labels_gp = gpar(fontsize = 15)
)

jpeg("TF_ranking.jpg", width=3000, height=1000, res=300)

library(circlize)
library(viridis)
ht <- Heatmap(as.matrix(t(toPlot[2:4])), name = "p.value",
        left_annotation = ha1, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c(viridis(100)[100],viridis(100)[75], viridis(100)[50], viridis(100)[25], viridis(100)[1] )),
        column_labels = toPlot$TF, heatmap_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 10), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE, cell_fun = cell_fun)

draw(ht, annotation_legend_list = list(dot_legend))

dev.off()

