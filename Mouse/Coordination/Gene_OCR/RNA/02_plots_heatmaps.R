###########################################
####     ATAC-Seq pipeline - DOUBLE    ####
####            04_heatmaps.R          ####
####              HUMAN                ####
###########################################

# Genes to plot 
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")

toPlot <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# Normalized matrix 
# regressing-out gsva for plotting 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
log_cpm_noBatch <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)

# Load fish data 

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
load("voom_to_dea.RData")
fish_rna <- fish_rna

#heatmap annotation
library(ComplexHeatmap)
library(RColorBrewer)

colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))
fish_anno <- fish_rna


fish_anno$translocation <- as.character(fish_anno$trans_4_14)
fish_anno[fish_anno$trans_14_16 == "t(14;16)",]$translocation <- "t(14;16)"
fish_anno$translocation <- as.character(fish_anno$translocation)
fish_anno[fish_anno$translocation == '',]$translocation <- "neutral"

fish_anno$qall <- as.character(fish_anno$qall)
fish_anno[fish_anno$qall == '',]$qall <- "neutral"
fish_anno$del1p <- as.character(fish_anno$del1p)
fish_anno[fish_anno$del1p == '',]$del1p <- "neutral"
fish_anno[fish_anno$del1p == "1p del",]$del1p <- "del(1p)"
fish_anno$del17p <- as.character(fish_anno$del17p)
fish_anno[fish_anno$del17p == '',]$del17p <- "neutral"
fish_anno[fish_anno$del17p == "17p del",]$del17p <- "del(17p)"
fish_anno$qall <- gsub("1q aber", "1q aberration", fish_anno$qall)

color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("HC", "MGUS", "SMM", "MM")

trans <- c("#ED90A4", "#C699E7","white")
names(trans) <- c("t(4;14)", "t(14;16)", "neutral")

amp1q <- c("#C0AB52", "white")
names(amp1q) <- c("1q aberration", "neutral")

del17 <- c("#4FBF85", "white")
names(del17) <- c("del(17p)", "neutral")

del1p <- c("#28BBD7", "white")
names(del1p) <- c("del(1p)", "neutral")

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
dim(count_table2) #1298


data_to_plot2 <- t(scale(t(count_table2)))
dim(data_to_plot2)



library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")
pdf("rna_dea_heatmap.pdf", width=10, height=10)
Heatmap(data_to_plot2, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()

# Plot heatmap with the clusters 

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

MGUS_HC_up <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==1,])]
SMM_HC_up <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==1,])]
MM_HC_up <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])]

MGUS_HC_down <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==-1,])]
SMM_HC_down <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==-1,])]
MM_HC_down <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])]

# UP in all
all_up <- intersect(MGUS_HC_up, intersect(SMM_HC_up, MM_HC_up))

# Down in all
all_down <- intersect(MGUS_HC_down, intersect(SMM_HC_down, MM_HC_down))

# Up in SMM and MM
up_SMM_MM <- intersect(SMM_HC_up, MM_HC_up)
up_SMM_MM <- up_SMM_MM[!up_SMM_MM %in% all_up]

# Down in SMM and MM
down_SMM_MM <- intersect(SMM_HC_down, MM_HC_down)
down_SMM_MM <- down_SMM_MM[!down_SMM_MM %in% all_down]

# Unique UP MM
unique_MM_up <- MM_HC_up[!MM_HC_up%in%c(MGUS_HC_up, SMM_HC_up)]

# Unique DOWN MM
unique_MM_down <- MM_HC_down[!MM_HC_down%in%c(MGUS_HC_down, SMM_HC_down)]

rest <- setdiff(c(MGUS_HC_up, SMM_HC_up, MM_HC_up, MGUS_HC_down, SMM_HC_down, MM_HC_down), c(all_up, all_down, up_SMM_MM, down_SMM_MM, unique_MM_up, unique_MM_down))

# Save 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA")
saveRDS(all_up, "all_up.RDS")
saveRDS(all_down, "all_down.RDS")
saveRDS(up_SMM_MM, "up_SMM_MM.RDS")
saveRDS(down_SMM_MM, "down_SMM_MM.RDS")
saveRDS(unique_MM_up, "unique_MM_up.RDS")
saveRDS(unique_MM_down, "unique_MM_down.RDS")
saveRDS(rest, "rest.RDS")


# Add cluster information to the heatmap
count_table2$cluster <- "rest"
count_table2$cluster[rownames(count_table2) %in% all_up] <- "constant up"
count_table2$cluster[rownames(count_table2) %in% all_down] <- "constant down"
count_table2$cluster[rownames(count_table2) %in% up_SMM_MM] <- "SMM & MM up"
count_table2$cluster[rownames(count_table2) %in% down_SMM_MM] <- "SMM & MM down"
count_table2$cluster[rownames(count_table2) %in% unique_MM_up] <- "MM up"
count_table2$cluster[rownames(count_table2) %in% unique_MM_down] <- "MM down"

data_to_plot2 <- t(scale(t(count_table2[,1:(ncol(count_table2)-1)])))

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")
pdf("rna_dea_heatmap_clusters.pdf", width=16, height=14)
Heatmap(data_to_plot2, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, row_split = factor(count_table2$cluster, levels = c( "constant down", "SMM & MM down", "MM down",  "MM up", "SMM & MM up", "constant up", "rest")),
        column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0.5, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()