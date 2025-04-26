###########################################
####     ATAC-Seq pipeline - DOUBLE    ####
####            04_heatmaps.R          ####
####              HUMAN                ####
###########################################

# Genes to plot 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
dea_results <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
CyclinD1_MGUS_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS_HC <- rownames(dea_results[dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS_HC <- rownames(dea_results[dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS_HC <- rownames(dea_results[dea_results$sig_MIc_MGUS != 0,])
CyclinD1_MM_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM_HC <- rownames(dea_results[dea_results$sig_Mmset_MM != 0,])
Trp53_MM_HC <- rownames(dea_results[dea_results$sig_Trp53_MM != 0,])
MIc_MM_HC <- rownames(dea_results[dea_results$sig_MIc_MM != 0,])

toPlot <- unique(c(CyclinD1_MGUS_HC, Mmset_MGUS_HC, Trp53_MGUS_HC, MIc_MGUS_HC,
                            CyclinD1_MM_HC, Mmset_MM_HC, Trp53_MM_HC, MIc_MM_HC))


# Normalized matrix 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")
log_cpm <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)
colnames(log_cpm) <- gsub("X", "", colnames(log_cpm))


#metadata
metadata <- read.table("Col_Data.csv",sep=",", header=TRUE, check.names=FALSE) 
metadata <- metadata[metadata$Sample_ID %in% colnames(log_cpm),]
log_cpm <- log_cpm[,metadata$Sample_ID]


#heatmap annotation
library(ComplexHeatmap)
library(RColorBrewer)


# Labels colors
color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")


color_models <- c("#368716", "#CAB2DC",  "#FB9A99", "#CAF2B0", "#1F78B4")
names(color_models) <-c ("Control","CyclinD1_BIc", "Mmset_BIc", "Trp53_BIc", "MIc")



# library(circlize)
# col_fun = colorRamp2(c(min(fish_anno$pc_infilt), max(fish_anno$pc_infilt)), c("blue","red"))

metadata$Model <- factor(metadata$Model, levels = c("Control","CyclinD1_BIc", "Mmset_BIc", "Trp53_BIc", "MIc"))

ha1 <- HeatmapAnnotation(
  Stage = metadata$Stage,
  Model = metadata$Model,
  col = list(Stage = as.factor(color_stage),
        Model = as.factor(color_models)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- log_cpm[toPlot,]
dim(count_table2)


data_to_plot2 <- t(scale(t(count_table2)))
dim(data_to_plot2)



library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("rna_dea_heatmap.pdf", width=10, height=8)
# order samples by model
Heatmap(data_to_plot2, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(metadata$Combined, levels = c("Control", "CyclinD1_BIc_MGUS", "Mmset_BIc_MGUS", "Trp53_BIc_MGUS", "MIc_MGUS",
                                                                                  "CyclinD1_BIc_MM", "Mmset_BIc_MM", "Trp53_BIc_MM", "MIc_MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()