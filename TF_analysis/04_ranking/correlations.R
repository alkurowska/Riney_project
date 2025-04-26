# correlation
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
TF_data <- read.table("TFs.txt", header = TRUE)
head(TF_data)

library(ggrepel)
library(ggplot2)

colnames(TF_data)
toPlot <- data.frame(TF_data$OR_MGUS_HC, TF_data$OR_SMM_HC, TF_data$OR_MM_HC,
                    TF_data$OR_inv_MGUS_HC, TF_data$OR_inv_SMM_HC, TF_data$OR_inv_MM_HC,
                    TF_data$DE_logFC_MGUS_HC, TF_data$DE_logFC_SMM_HC, TF_data$DE_logFC_MM_HC, 
                    TF_data$chromVAR_logFC_MGUS_HC, TF_data$chromVAR_logFC_SMM_HC, TF_data$chromVAR_logFC_MM_HC,
                    TF_data$abs_NES_MGUS_HC, TF_data$abs_NES_SMM_HC, TF_data$abs_NES_MM_HC)
colnames(toPlot) <- c("MGUS_HC_OR", "SMM_HC_OR", "MM_HC_OR", 
                      "MGUS_HC_OR_inv", "SMM_HC_OR_inv", "MM_HC_OR_inv",
                      "MGUS_HC_DE_logFC", "SMM_HC_DE_logFC", "MM_HC_DE_logFC",
                      "MGUS_HC_chromVAR_logFC", "SMM_HC_chromVAR_logFC", "MM_HC_chromVAR_logFC",
                      "MGUS_HC_NES_abs", "SMM_HC_NES_abs", "MM_HC_NES_abs")

# Add rownames, if duplicated add ".number"
rownames(toPlot) <- make.unique(TF_data$TFs)


Transition <- sapply(strsplit(colnames(toPlot),"_"), `[`, 1)
Test <- gsub("MGUS_HC_", "", colnames(toPlot)) 
Test <- gsub("SMM_HC_", "", Test)
Test <- gsub("MM_HC_", "", Test)

metadata <- data.frame(Transition, Test)
metadata$Test <- factor(metadata$Test, levels = c("OR", 
                        "OR_inv",
                        "DE_logFC",  
                        "chromVAR_logFC",
                        "NES_abs"))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

color_trans <- c("#D27D46", "#D92F5E", "#8D4E85") 
names(color_trans) <- c("MGUS", "SMM", "MM")


color_test <- c("bisque3", "bisque4", "darkslategray", 
                "azure4", #"darkslategray1", "darkturquoise", 
                "deepskyblue3")
names(color_test) <- unique(Test)


ha1 <- rowAnnotation(
  Transition = metadata$Transition,
  Test = metadata$Test,
  col = list(Transition = as.factor(color_trans),
             Test = as.factor(color_test)),
  simple_anno_size = unit(1, "cm"),
  na_col = "#D9D9D9", annotation_name_gp= gpar(fontsize = 20),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
  show_annotation_name = TRUE)

ha <- HeatmapAnnotation(
  Transition = metadata$Transition,
  Test = metadata$Test,
  col = list(Transition = as.factor(color_trans),
             Test = as.factor(color_test)),
  simple_anno_size = unit(1, "cm"),
  na_col = "#D9D9D9", annotation_name_gp= gpar(fontsize = 20),
  show_legend = FALSE,
  annotation_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
  show_annotation_name = TRUE)

getwd()
jpeg("TF_heatmap_trans.jpg", width=8000, height=3000, res=300)

library(circlize)

Heatmap(as.matrix(t(scale(toPlot))), name = "z-score",
        left_annotation = ha1, col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        column_labels = rownames(toPlot), heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE)

dev.off()


jpeg("TF_heatmap_test.jpg", width=8000, height=3000, res=300)

Heatmap(as.matrix(t(scale(toPlot))), name = "z-score",
        left_annotation = ha1, col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        column_labels = rownames(toPlot), heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Test, levels = c("OR", 
                        "OR_inv",
                        "DE_logFC",  
                        "chromVAR_logFC",
                        "NES_abs")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE)

dev.off()



# correlation plot 
# create a data frame with correlations between all columns of toPlot 
correlations <- cor(toPlot, method = "spearman")
# plot a heatmap of the correlations

jpeg("correlation_heatmap_trans.jpg", width=7000, height=6000, res=300)   
Heatmap(correlations, name = "Correlation", left_annotation = ha1, top_annotation = ha, col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontsize = 15), 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = TRUE, row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0, column_title_rot = 0, cluster_row_slices = F, cluster_column_slices = F,
        column_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        cluster_columns = F, cluster_rows = F, use_raster=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 15))
})
dev.off()

jpeg("correlation_heatmap_test.jpg", width=7000, height=6000, res=300)   
Heatmap(correlations, name = "Correlation", left_annotation = ha1, top_annotation = ha, col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontsize = 15), 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = TRUE, row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        row_split = factor(metadata$Test, levels = c("OR", "OR_inv", "DE_logFC", "chromVAR_logFC", "NES_abs")),
        row_title_rot = 0, column_title_rot = 0, cluster_row_slices = F, cluster_column_slices = F,
        column_split = factor(metadata$Test, levels = c("OR", "OR_inv", "DE_logFC", "chromVAR_logFC", "NES_abs")),
        cluster_columns = F, cluster_rows = F, use_raster=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 15))
})
dev.off()



# Plot with NES pos and NES neg
toPlot <- data.frame(TF_data$OR_MGUS_HC, TF_data$OR_SMM_HC, TF_data$OR_MM_HC,
                    TF_data$OR_inv_MGUS_HC, TF_data$OR_inv_SMM_HC, TF_data$OR_inv_MM_HC,
                    TF_data$DE_logFC_MGUS_HC, TF_data$DE_logFC_SMM_HC, TF_data$DE_logFC_MM_HC, 
                    TF_data$chromVAR_logFC_MGUS_HC, TF_data$chromVAR_logFC_SMM_HC, TF_data$chromVAR_logFC_MM_HC,
                    TF_data$pos_NES_MGUS_HC, TF_data$pos_NES_SMM_HC, TF_data$pos_NES_MM_HC,
                    TF_data$neg_NES_MGUS_HC, TF_data$neg_NES_SMM_HC, TF_data$neg_NES_MM_HC,
                                        TF_data$abs_NES_MGUS_HC, TF_data$abs_NES_SMM_HC, TF_data$abs_NES_MM_HC)
colnames(toPlot) <- c("MGUS_HC_OR", "SMM_HC_OR", "MM_HC_OR", 
                      "MGUS_HC_OR_inv", "SMM_HC_OR_inv", "MM_HC_OR_inv",
                      "MGUS_HC_DE_logFC", "SMM_HC_DE_logFC", "MM_HC_DE_logFC",
                      "MGUS_HC_chromVAR_logFC", "SMM_HC_chromVAR_logFC", "MM_HC_chromVAR_logFC",
                      "MGUS_HC_NES_pos", "SMM_HC_NES_pos", "MM_HC_NES_pos",
                      "MGUS_HC_NES_neg", "SMM_HC_NES_neg", "MM_HC_NES_neg",
                      "MGUS_HC_NES_abs", "SMM_HC_NES_abs", "MM_HC_NES_abs")

# Add rownames, if duplicated add ".number"
rownames(toPlot) <- make.unique(TF_data$TFs)


Transition <- sapply(strsplit(colnames(toPlot),"_"), `[`, 1)
Test <- gsub("MGUS_HC_", "", colnames(toPlot)) 
Test <- gsub("SMM_HC_", "", Test)
Test <- gsub("MM_HC_", "", Test)

metadata <- data.frame(Transition, Test)
metadata$Test <- factor(metadata$Test, levels = c("OR", 
                        "OR_inv",
                        "DE_logFC",  
                        "chromVAR_logFC",
                        "NES_pos",
                        "NES_neg",
                        "NES_abs"))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

color_trans <- c("#D27D46", "#D92F5E", "#8D4E85") 
names(color_trans) <- c("MGUS", "SMM", "MM")


color_test <- c("bisque3", "bisque4", "darkslategray", 
                "azure4", "darkslategray1", "darkturquoise", 
                "deepskyblue3")
names(color_test) <- unique(Test)


ha1 <- rowAnnotation(
  Transition = metadata$Transition,
  Test = metadata$Test,
  col = list(Transition = as.factor(color_trans),
             Test = as.factor(color_test)),
  simple_anno_size = unit(1, "cm"),
  na_col = "#D9D9D9", annotation_name_gp= gpar(fontsize = 20),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
  show_annotation_name = TRUE)

ha <- HeatmapAnnotation(
  Transition = metadata$Transition,
  Test = metadata$Test,
  col = list(Transition = as.factor(color_trans),
             Test = as.factor(color_test)),
  simple_anno_size = unit(1, "cm"),
  na_col = "#D9D9D9", annotation_name_gp= gpar(fontsize = 20),
  show_legend = FALSE,
  annotation_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
  show_annotation_name = TRUE)


jpeg("TF_heatmap_trans_all.jpg", width=8000, height=3000, res=300)

library(circlize)

Heatmap(as.matrix(t(scale(toPlot))), name = "z-score",
        left_annotation = ha1, col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        column_labels = rownames(toPlot), heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE)

dev.off()


jpeg("TF_heatmap_test_all.jpg", width=8000, height=3000, res=300)

Heatmap(as.matrix(t(scale(toPlot))), name = "z-score",
        left_annotation = ha1, col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        column_labels = rownames(toPlot), heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 4), 
        cluster_columns = TRUE, 
        row_split = factor(metadata$Test, levels = c("OR", 
                        "OR_inv",
                        "DE_logFC",  
                        "chromVAR_logFC",
                        "NES_pos",
                        "NES_neg",
                        "NES_abs")),
        row_title_rot = 0,
        cluster_row_slices = FALSE, row_gap = unit(2, "mm"),
        cluster_rows = F, use_raster=TRUE)

dev.off()



# correlation plot 
# create a data frame with correlations between all columns of toPlot 
# remove rows with NA
toPlot <- toPlot[complete.cases(toPlot), ]
correlations <- cor(toPlot, method = "spearman")
# plot a heatmap of the correlations

jpeg("correlation_heatmap_trans_all.jpg", width=7000, height=6000, res=300)   
Heatmap(correlations, name = "Correlation", left_annotation = ha1, top_annotation = ha, col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontsize = 15), 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = TRUE, row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        row_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        row_title_rot = 0, column_title_rot = 0, cluster_row_slices = F, cluster_column_slices = F,
        column_split = factor(metadata$Transition, levels = c("MGUS", "SMM", "MM")),
        cluster_columns = F, cluster_rows = F, use_raster=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 15))
})
dev.off()

jpeg("correlation_heatmap_test_all.jpg", width=7000, height=6000, res=300)   
Heatmap(correlations, name = "Correlation", left_annotation = ha1, top_annotation = ha, col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        column_names_gp = gpar(fontsize = 15), row_names_gp = gpar(fontsize = 15), 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20, fontface = "bold")),
        show_column_names = TRUE, show_row_names = TRUE, row_gap = unit(3, "mm"), column_gap = unit(3, "mm"),
        row_split = factor(metadata$Test, levels = c("OR", "OR_inv", "DE_logFC",  "chromVAR_logFC", "NES_pos", "NES_neg", "NES_abs")),
        row_title_rot = 0, column_title_rot = 0, cluster_row_slices = F, cluster_column_slices = F,
        column_split = factor(metadata$Test, levels = c("OR", "OR_inv", "DE_logFC", "chromVAR_logFC", "NES_pos", "NES_neg", "NES_abs")),
        cluster_columns = F, cluster_rows = F, use_raster=TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", correlations[i, j]), x, y, gp = gpar(fontsize = 15))
})
dev.off()
