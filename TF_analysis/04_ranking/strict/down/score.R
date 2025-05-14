##### READ DATA #####
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking/strict/down")
MGUS <- read.table("Rank_MGUS.txt",header = T, sep = "\t")
SMM <- read.table("Rank_SMM.txt",header = T, sep = "\t")
MM <- read.table("Rank_MM.txt",header = T, sep = "\t")

# sort by name
MGUS <- MGUS[order(MGUS[,1]),]
SMM <- SMM[order(SMM[,1]),]
MM <- MM[order(MM[,1]),]
summary(rownames(MGUS) == rownames(SMM)) # TRUE
summary(rownames(MGUS) == rownames(MM)) # TRUE


# FDR 
MGUS$p.adjust <- p.adjust(MGUS$Score, method = "fdr")
SMM$p.adjust <- p.adjust(SMM$Score, method = "fdr")
MM$p.adjust <- p.adjust(MM$Score, method = "fdr")


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
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
Pioneer <- as.character( read.table("Pioneers.txt",header = F)[,1])
Crisp <- as.character(read.table("Crispr.txt",header = F)[,1])

# Are any of the essential TFs in the list?
Pioneer[(Pioneer %in% toPlot$TF)]
Crisp[(Crisp %in% toPlot$TF)] #"MEF2D" "MEF2C"


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

# Load expression data

# Normalized matrix 
# regressing-out gsva for plotting 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
log_cpm_noBatch <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)
colnames(log_cpm_noBatch) <- gsub("X", "", colnames(log_cpm_noBatch))

# Load fish data 

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
load("voom_to_dea.RData")
fish_rna <- fish_rna
gene_anno <- gene_anno

TFs <- c(toPlot$TF)
# remove . and a number after the dot in TFs names
TFs <- gsub("\\..*", "", TFs)
TFs <- as.data.frame(TFs)
TFs$gene_id <- gene_anno$gene_id[match(TFs$TFs, gene_anno$gene_name)]

TFs$avg_expr_MGUS <- rowMeans(log_cpm_noBatch[TFs$gene_id, fish_rna$`Sample RNAseq`[fish_rna$Stage == "MGUS"]], na.rm = TRUE)
TFs$avg_expr_SMM  <- rowMeans(log_cpm_noBatch[TFs$gene_id, fish_rna$`Sample RNAseq`[fish_rna$Stage == "SMM"]], na.rm = TRUE)
TFs$avg_expr_MM   <- rowMeans(log_cpm_noBatch[TFs$gene_id, fish_rna$`Sample RNAseq`[fish_rna$Stage == "MM"]], na.rm = TRUE)

# Combine into a data frame for annotation
avg_expr_df <- data.frame(
  MGUS = TFs$avg_expr_MGUS,
  SMM  = TFs$avg_expr_SMM,
  MM   = TFs$avg_expr_MM
)

top_anno <- HeatmapAnnotation(
  "AvgExpr_MGUS" = anno_barplot(avg_expr_df$MGUS, gp = gpar(fill = color_trans["MGUS"]), border = FALSE, bar_width = 0.6,  axis = TRUE, axis_param = list(gp = gpar(fontsize = 4))),  # Set y-axis label font size here,
  "AvgExpr_SMM"  = anno_barplot(avg_expr_df$SMM,  gp = gpar(fill = color_trans["SMM"]),  border = FALSE, bar_width = 0.6,  axis = TRUE, axis_param = list(gp = gpar(fontsize = 4))),  
  "AvgExpr_MM"   = anno_barplot(avg_expr_df$MM,   gp = gpar(fill = color_trans["MM"]),   border = FALSE, bar_width = 0.6,  axis = TRUE, axis_param = list(gp = gpar(fontsize = 4))),  
  annotation_name_side = "left", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10)
)

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking/strict/down")
jpeg("TF_ranking.jpg", width=3000, height=1500, res=300)

library(viridis)
library(circlize)


ht <- Heatmap(as.matrix(t(toPlot[2:4])), name = "p.value", top_annotation = top_anno, 
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



# add a logFC instead 
# load dea results 

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_results_noscore <- read.table("rna_dea_results_noscore.txt", header = TRUE, sep = "\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header = TRUE, sep = "\t")


TFs$logFC_MGUS_HC <- dea_results[TFs$gene_id, ]$logFC_MGUS_HC
TFs$logFC_SMM_HC  <- dea_results[TFs$gene_id, ]$logFC_SMM_HC
TFs$logFC_MM_HC   <- dea_results_noscore[TFs$gene_id, ]$logFC_MM_HC

# Combine into a data frame for annotation
logFC <- data.frame(
  MGUS = TFs$logFC_MGUS_HC,
  SMM  = TFs$logFC_SMM_HC,
  MM   = TFs$logFC_MM_HC
)
mat_logFC <- as.matrix(logFC)

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking/strict/down")
jpeg("TF_ranking_logFC.jpg", width=3000, height=1000, res=300)

library(viridis)
library(circlize)


ht <- Heatmap(as.matrix(t(toPlot[2:4])), name = "p.value",
        left_annotation = ha1, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c(viridis(100)[100],viridis(100)[75], viridis(100)[50], viridis(100)[25], viridis(100)[1] )),
        column_labels = toPlot$TF, heatmap_legend_param = list(labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 15, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 10), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE, cell_fun = cell_fun)


# Get the column order after clustering
ht_drawn <- draw(ht, heatmap_legend_side = "right", annotation_legend_list = list(dot_legend), merge_legend = TRUE, show_heatmap_legend = FALSE)
col_order <- column_order(ht_drawn)

top_heatmap <- Heatmap(
  t(mat_logFC),
  name = "logFC",
  col = colorRamp2(
    c(-5, 0, 5),
     c("#4475B3", "white", "#D7342A")
  ),
  show_row_names = TRUE, 
  # show row names on the left side
  row_names_side = "left",
  # row labels define
  row_labels = c("MGUS vs HC", "SMM vs HC", "MM vs HC"),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,  # <--- turn off clustering
  column_order = col_order,
  height = unit(1.2, "cm"),
  column_title = NULL,
  row_title = NULL,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 4),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)


# Combine the top heatmap and the main heatmap
ht_list <- top_heatmap %v% ht

draw(ht_list, annotation_legend_list = list(dot_legend))


dev.off()


TFs