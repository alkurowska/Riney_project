###########################################
####     ATAC-Seq pipeline - DOUBLE    ####
####            04_heatmaps.R          ####
####              HUMAN                ####
###########################################

# Genes to plot 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")
toPlot <- unique(c(MGUS_HC_filter, SMM_HC_filter, MM_HC_filter))

# Normalized matrix 
# regressing-out gsva for plotting 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
log_cpm_noBatch <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)

# Load fish data 

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
load("voom_to_dea.RData")
fish_atac <- fish_atac

#heatmap annotation
library(ComplexHeatmap)
library(RColorBrewer)

colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))
fish_anno <- fish_atac


fish_anno$translocation <- as.character(fish_anno$trans_4_14)
fish_anno[fish_anno$trans_14_16 == "t(14;16)",]$translocation <- "t(14;16)"
fish_anno$translocation <- as.character(fish_anno$translocation)
fish_anno[fish_anno$translocation == '',]$translocation <- "neutral"

fish_anno$qall <- as.character(fish_anno$qall)
fish_anno[fish_anno$qall == '',]$qall <- "neutral"
fish_anno$del1p <- as.character(fish_anno$del1p)
fish_anno[fish_anno$del1p == '',]$del1p <- "neutral"
fish_anno$del17p <- as.character(fish_anno$del17p)
fish_anno[fish_anno$del17p == '',]$del17p <- "neutral"
fish_anno$qall <- gsub("1q aber", "1q aberration", fish_anno$qall)

color_stage <- c("#33A02C", "#FFBF00","#FF7F00", "#DC143C") 
names(color_stage) <- c("HC", "MGUS", "SMM", "MM")

trans <- c("#ED90A4", "#C699E7","white")
names(trans) <- c("t(4;14)", "t(14;16)", "neutral")

amp1q <- c("#C0AB52", "white")
names(amp1q) <- c("1q aberration", "neutral")

del17 <- c("#4FBF85", "white")
names(del17) <- c("17p del", "neutral")

del1p <- c("#28BBD7", "white")
names(del1p) <- c("1p del", "neutral")

# library(circlize)
# col_fun = colorRamp2(c(min(fish_anno$pc_infilt), max(fish_anno$pc_infilt)), c("blue","red"))


ha1 <- HeatmapAnnotation(
  Stage = fish_anno$Stage,
  `1q aberration` = fish_anno$qall,
   Translocations = fish_anno$translocation,
  `17p deletion` = fish_anno$del17p,
  `1p deletion` = fish_anno$del1p,
  col = list(Stage = as.factor(color_stage),
         `1q aberration` = as.factor(amp1q),
             Translocations = as.factor(trans),
             `17p deletion` = as.factor(del17),
             `1p deletion` = as.factor(del1p)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- log_cpm_noBatch[toPlot,]
dim(count_table2)


data_to_plot2 <- t(scale(t(count_table2)))
dim(data_to_plot2)



library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/PLOTS")
pdf("atac_dea_heatmap.pdf", width=10, height=12)
Heatmap(data_to_plot2, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()