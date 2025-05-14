###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

library(UpSetR)

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

CyclinD1_MGUS_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS_HC <- rownames(dea_results[dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS_HC <- rownames(dea_results[dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS_HC <- rownames(dea_results[dea_results$sig_MIc_MGUS != 0,])
CyclinD1_MM_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM_HC <- rownames(dea_results[dea_results$sig_Mmset_MM != 0,])
Trp53_MM_HC <- rownames(dea_results[dea_results$sig_Trp53_MM != 0,])
MIc_MM_HC <- rownames(dea_results[dea_results$sig_MIc_MM != 0,])

# plot without dividing per contrast
CyclinD1 <- unique(c(CyclinD1_MGUS_HC, CyclinD1_MM_HC))
Mmset <- unique(c(Mmset_MGUS_HC, Mmset_MM_HC))
Trp53 <- unique(c(Trp53_MGUS_HC, Trp53_MM_HC))
MIc <- unique(c(MIc_MGUS_HC, MIc_MM_HC))


# LOAD FINAL human genes
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")
final_genes <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# LOAD all of the human genes
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
load("norm_to_voom.RData")
gene_anno_human <- gene_anno
gene_anno_human <- gene_anno_human[gene_anno_human$gene_id%in%rownames(d1_norm$counts),]
dim(gene_anno_human) # 24256 genes

all_genes <- gene_anno_human$gene_id

# Load orthologs 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
orthologs <- read.table("human_mouse_orthologs.txt", sep = "\t", header = TRUE)

all_genes <- unique(orthologs$Mouse_ensembl_gene_id[orthologs$Human_ensembl_gene_id %in% all_genes])

Human_ocr_genes <- unique(orthologs$Mouse_ensembl_gene_id[orthologs$Human_ensembl_gene_id %in% final_genes])

sets <- unique(c(CyclinD1, Mmset, Trp53, MIc, #all_genes, 
Human_ocr_genes))
# Create a binary matrix of membership

genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = genes,
  CyclinD1 = genes %in% CyclinD1,
  Mmset = genes %in% Mmset,
  Trp53 = genes %in% Trp53,
  MIc = genes %in% MIc,
  #Human_all_genes = genes %in% all_genes,
  OCR_gene = genes %in% Human_ocr_genes
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("CyclinD1", "Mmset", "Trp53", "MIc", #"All orthologs", 
"OCR-gene orthologs")


# Generate the UpSet plot

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_general_no_contrast_orthologs.png", width=3000, height=1500, res=300)
upset(upset_data, sets = c("CyclinD1", "Mmset", "Trp53", "MIc", #"All orthologs", 
"OCR-gene orthologs"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

# find orthologs overlaping with any of the model
all_mouse <- unique(c(CyclinD1, Mmset, Trp53, MIc))
orthologs_any <- intersect(Human_ocr_genes, all_mouse)
orthologs_any <- orthologs[orthologs$Mouse_ensembl_gene_id %in% orthologs_any,]
write.table(orthologs_any, file = "orthologs_any.txt", sep = "\t", quote = FALSE, row.names = FALSE)


genes85 <- intersect(intersect(intersect(CyclinD1, Mmset), intersect(Trp53, MIc)), Human_ocr_genes)
orthologs_85 <- orthologs[orthologs$Mouse_ensembl_gene_id %in% genes85,]
write.table(orthologs_85, file = "orthologs_85.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dim(orthologs_85) # 85 genes
head(orthologs_85)

colnames(dea_results) # 4 columns
# Add differential results
toPlot_85 <- merge(orthologs_85, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA/")
human_dea <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")
toPlot_85 <- merge(toPlot_85, human_dea[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
human_dea_no_score <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")
toPlot_85$sig_MM_HC <- human_dea_no_score[toPlot_85$Human_ensembl_gene_id,]$sig_MM_HC

colnames(toPlot_85) <- gsub("sig_", "", colnames(toPlot_85))
library(stringr)
metadata <- as.data.frame(str_split_fixed(colnames(toPlot_85)[3:ncol(toPlot_85)], "_", n = 2))
colnames(metadata) <- c("Model", "Contrast")
metadata$Model <- gsub("MGUS", "Human", metadata$Model)
metadata$Model <- gsub("SMM", "Human", metadata$Model)
metadata$Model <- gsub("MM", "Human", metadata$Model)
metadata$Contrast[9:11] <- c("MGUS", "SMM", "MM")

rownames(metadata) <- colnames(toPlot_85)[3:ncol(toPlot_85)]


color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")

color_models <- c("#368716", "#CAB2DC",  "#FB9A99", "#CAF2B0", "#1F78B4","yellow2")
names(color_models) <-c ("Control","CyclinD1", "Mmset", "Trp53", "MIc", "Human")


metadata$Model <- factor(metadata$Model, levels = c("CyclinD1", "Mmset", "Trp53", "MIc", "Human"))

library(ComplexHeatmap)
library(circlize)

ha1 <- HeatmapAnnotation(
  Contrast = metadata$Contrast,
  Model = metadata$Model,
  col = list(Contrast = as.factor(color_stage),
        Model = as.factor(color_models)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- toPlot_85[,3:ncol(toPlot_85)]
dim(count_table2)


library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs85_heatmap.pdf", width=5, height=8)
# order samples by model
Heatmap(as.matrix(count_table2), name = "Differential Expression", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(count_table2)*unit(5, "mm"), #height = nrow(count_table2)*unit(5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()

# Plot logFC instead
# Add differential results
toPlot_85 <- merge(orthologs_85, dea_results[,c("logFC_CyclinD1_MGUS", "logFC_CyclinD1_MM", "logFC_MIc_MM", "logFC_MIc_MGUS", "logFC_Mmset_MM", "logFC_Mmset_MGUS","logFC_Trp53_MM", "logFC_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
toPlot_85 <- merge(toPlot_85, human_dea[,c("logFC_MGUS_HC", "logFC_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
toPlot_85$logFC_MM_HC <- human_dea_no_score[toPlot_85$Human_ensembl_gene_id,]$logFC_MM_HC

colnames(toPlot_85) <- gsub("logFC_", "", colnames(toPlot_85))
library(stringr)
metadata <- as.data.frame(str_split_fixed(colnames(toPlot_85)[3:ncol(toPlot_85)], "_", n = 2))
colnames(metadata) <- c("Model", "Contrast")
metadata$Model <- gsub("MGUS", "Human", metadata$Model)
metadata$Model <- gsub("SMM", "Human", metadata$Model)
metadata$Model <- gsub("MM", "Human", metadata$Model)
metadata$Contrast[9:11] <- c("MGUS", "SMM", "MM")

rownames(metadata) <- colnames(toPlot_85)[3:ncol(toPlot_85)]


color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")

color_models <- c("#368716", "#CAB2DC",  "#FB9A99", "#CAF2B0", "#1F78B4","yellow2")
names(color_models) <-c ("Control","CyclinD1", "Mmset", "Trp53", "MIc", "Human")


metadata$Model <- factor(metadata$Model, levels = c("CyclinD1", "Mmset", "Trp53", "MIc", "Human"))

library(ComplexHeatmap)
library(circlize)

ha1 <- HeatmapAnnotation(
  Contrast = metadata$Contrast,
  Model = metadata$Model,
  col = list(Contrast = as.factor(color_stage),
        Model = as.factor(color_models)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- toPlot_85[,3:ncol(toPlot_85)]
dim(count_table2)

data_to_plot <- t(scale(t(count_table2)))

library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs85_heatmap_logFC.pdf", width=5, height=8)
# order samples by model
Heatmap(data_to_plot, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(data_to_plot)*unit(5, "mm"), 
        #height = nrow(data_to_plot)*unit(5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()