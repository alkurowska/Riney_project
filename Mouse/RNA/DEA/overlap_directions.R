###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

# Load the HUMAN DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
human_dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
human_dea_results_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load overlaping orthologs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
orthologs_85 <- read.table("orthologs_85.txt", sep = "\t", header = TRUE)

# ADD dynamics 
# Add differential results
orthologs_85 <- merge(orthologs_85, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
orthologs_85 <- merge(orthologs_85, human_dea_results[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
orthologs_85$sig_MM_HC <- human_dea_results_noscore[orthologs_85$Human_ensembl_gene_id,]$sig_MM_HC

# Keep rows only with 0 and 1
orthologs_85_up <- orthologs_85[orthologs_85$sig_CyclinD1_MGUS %in% c(0, 1) & orthologs_85$sig_CyclinD1_MM %in% c(0, 1) & orthologs_85$sig_MIc_MM %in% c(0, 1) & orthologs_85$sig_MIc_MGUS %in% c(0, 1) & orthologs_85$sig_Mmset_MM %in% c(0, 1) & orthologs_85$sig_Mmset_MGUS %in% c(0, 1) & orthologs_85$sig_Trp53_MM %in% c(0, 1) & orthologs_85$sig_Trp53_MGUS %in% c(0, 1) & orthologs_85$sig_MGUS_HC %in% c(0, 1) & orthologs_85$sig_SMM_HC %in% c(0, 1) & orthologs_85$sig_MM_HC %in% c(0, 1),]
# remove rows with 0 for all human contrasts
orthologs_85_up <- orthologs_85_up[!(orthologs_85_up$sig_MGUS_HC == 0 & orthologs_85_up$sig_SMM_HC == 0 & orthologs_85_up$sig_MM_HC == 0),]

# Keep rows only with 0 and -1
orthologs_85_down <- orthologs_85[orthologs_85$sig_CyclinD1_MGUS %in% c(0, -1) & orthologs_85$sig_CyclinD1_MM %in% c(0, -1) & orthologs_85$sig_MIc_MM %in% c(0, -1) & orthologs_85$sig_MIc_MGUS %in% c(0, -1) & orthologs_85$sig_Mmset_MM %in% c(0, -1) & orthologs_85$sig_Mmset_MGUS %in% c(0, -1) & orthologs_85$sig_Trp53_MM %in% c(0, -1) & orthologs_85$sig_Trp53_MGUS %in% c(0, -1) & orthologs_85$sig_MGUS_HC %in% c(0, -1) & orthologs_85$sig_SMM_HC %in% c(0, -1) & orthologs_85$sig_MM_HC %in% c(0, -1),]
# remove rows with 0 for all human contrasts
orthologs_85_down <- orthologs_85_down[!(orthologs_85_down$sig_MGUS_HC == 0 & orthologs_85_down$sig_SMM_HC == 0 & orthologs_85_down$sig_MM_HC == 0),]


dim(orthologs_85_up) # 13
dim(orthologs_85_down) # 40
# Keep rows only with 0 and 1

# Add gene names 
setwd("/ibex/user/kurowsaa/RINEY/human/RNA_data")
gene_anno <- readRDS("Complete_GTF.RDS")

orthologs_85_up$gene_name <- gene_anno[match(orthologs_85_up$Human_ensembl_gene_id, gene_anno$gene_id),]$gene_name
orthologs_85_down$gene_name <- gene_anno[match(orthologs_85_down$Human_ensembl_gene_id, gene_anno$gene_id),]$gene_name

sort(unique(orthologs_85_up$gene_name))
sort(unique(orthologs_85_down$gene_name))


# List of essential TFs
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
Pioneer <- as.character( read.table("Pioneers.txt",header = F)[,1])
Crisp <- as.character(read.table("Crispr.txt",header = F)[,1])

Pioneer[(Pioneer %in% orthologs_85_up$gene_name)]
Crisp[(Crisp %in% orthologs_85_up$gene_name)]
Pioneer[(Pioneer %in% orthologs_85_down$gene_name)] #"CEBPA" "EBF1"
Crisp[(Crisp %in% orthologs_85_down$gene_name)] 

toPlot_up <- orthologs_85_up
colnames(toPlot_up) <- gsub("sig_", "", colnames(toPlot_up))
library(stringr)
metadata <- as.data.frame(str_split_fixed(colnames(toPlot_up)[3:(ncol(toPlot_up)-1)], "_", n = 2))
colnames(metadata) <- c("Model", "Contrast")
metadata$Model <- gsub("MGUS", "Human", metadata$Model)
metadata$Model <- gsub("SMM", "Human", metadata$Model)
metadata$Model <- gsub("MM", "Human", metadata$Model)
metadata$Contrast[9:11] <- c("MGUS", "SMM", "MM")

rownames(metadata) <- colnames(toPlot_up)[3:(ncol(toPlot_up)-1)]


color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")

color_models <- c("#368716", "#CAB2DC",  "#FB9A99", "#CAF2B0", "#1F78B4","yellow2")
names(color_models) <-c ("Control","CyclinD1", "Mmset", "Trp53", "MIc", "Human")


metadata$Model <- factor(metadata$Model, levels = c("CyclinD1", "Mmset", "Trp53", "MIc", "Human"))
metadata$Combined <- paste(metadata$Model, metadata$Contrast, sep = "_")
metadata$Combined <- gsub("Human_MGUS", "MGUS_HC", metadata$Combined)
metadata$Combined <- gsub("Human_SMM", "SMM_HC", metadata$Combined)
metadata$Combined <- gsub("Human_MM", "MM_HC", metadata$Combined)
metadata$Combined <- factor(metadata$Combined, levels = c( "MGUS_HC", "CyclinD1_MGUS", "Mmset_MGUS", "Trp53_MGUS", "MIc_MGUS","SMM_HC","CyclinD1_MM","Mmset_MM", "Trp53_MM", "MIc_MM", "MM_HC"))

library(ComplexHeatmap)
library(circlize)

ha1 <- HeatmapAnnotation(
  Contrast = metadata$Contrast,
  Model = metadata$Model,
  col = list(Contrast = as.factor(color_stage),
        Model = as.factor(color_models)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- toPlot_up[,3:(ncol(toPlot_up)-1)]
rownames(count_table2) <- toPlot_up$gene_name
dim(count_table2)
# Order columns to match metadata$Combined
column_ordered <- metadata$Combined[match(colnames(count_table2), rownames(metadata))]

library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs85_up_heatmap.pdf", width=6, height=3)
# order samples by model
Heatmap(as.matrix(count_table2), name = "Differential Expression", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), column_order = order(column_ordered),
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(count_table2)*unit(5, "mm"), #height = nrow(count_table2)*unit(5, "mm"),
        cluster_columns = FALSE, use_raster=TRUE)

dev.off()


toPlot_down <- orthologs_85_down
colnames(toPlot_down) <- gsub("sig_", "", colnames(toPlot_down))

count_table2 <- toPlot_down[,3:(ncol(toPlot_down)-1)]
rownames(count_table2) <- toPlot_down$gene_name
dim(count_table2)
# Order columns to match metadata$Combined
column_ordered <- metadata$Combined[match(colnames(count_table2), rownames(metadata))]


library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs85_down_heatmap.pdf", width=6, height=5)
# order samples by model
Heatmap(as.matrix(count_table2), name = "Differential Expression", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), column_order = order(column_ordered),
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(count_table2)*unit(5, "mm"), #height = nrow(count_table2)*unit(5, "mm"),
        cluster_columns = FALSE, use_raster=TRUE)

dev.off()




# Other genes 

# Load overlaping orthologs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
orthologs_any <- read.table("orthologs_any.txt", sep = "\t", header = TRUE)

# ADD dynamics 
# Add differential results
orthologs_any <- merge(orthologs_any, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
orthologs_any <- merge(orthologs_any, human_dea_results[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
orthologs_any$sig_MM_HC <- human_dea_results_noscore[orthologs_any$Human_ensembl_gene_id,]$sig_MM_HC

length(unique(orthologs_any$Mouse_ensembl_gene_id)) # 341
# remove rows from the list of 85
orthologs_any <- orthologs_any[!(orthologs_any$Mouse_ensembl_gene_id %in% orthologs_85$Mouse_ensembl_gene_id),]
length(unique(orthologs_any$Mouse_ensembl_gene_id)) # 256



# Keep rows only with 0 and 1
orthologs_any_up <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, 1) & orthologs_any$sig_CyclinD1_MM %in% c(0, 1) & orthologs_any$sig_MIc_MM %in% c(0, 1) & orthologs_any$sig_MIc_MGUS %in% c(0, 1) & orthologs_any$sig_Mmset_MM %in% c(0, 1) & orthologs_any$sig_Mmset_MGUS %in% c(0, 1) & orthologs_any$sig_Trp53_MM %in% c(0, 1) & orthologs_any$sig_Trp53_MGUS %in% c(0, 1) & orthologs_any$sig_MGUS_HC %in% c(0, 1) & orthologs_any$sig_SMM_HC %in% c(0, 1) & orthologs_any$sig_MM_HC %in% c(0, 1),]
#remove rows with 0 for all human contrasts
orthologs_any_up <- orthologs_any_up[!(orthologs_any_up$sig_MGUS_HC == 0 & orthologs_any_up$sig_SMM_HC == 0 & orthologs_any_up$sig_MM_HC == 0),]

# Keep rows only with 0 and -1
orthologs_any_down <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, -1) & orthologs_any$sig_CyclinD1_MM %in% c(0, -1) & orthologs_any$sig_MIc_MM %in% c(0, -1) & orthologs_any$sig_MIc_MGUS %in% c(0, -1) & orthologs_any$sig_Mmset_MM %in% c(0, -1) & orthologs_any$sig_Mmset_MGUS %in% c(0, -1) & orthologs_any$sig_Trp53_MM %in% c(0, -1) & orthologs_any$sig_Trp53_MGUS %in% c(0, -1) & orthologs_any$sig_MGUS_HC %in% c(0, -1) & orthologs_any$sig_SMM_HC %in% c(0, -1) & orthologs_any$sig_MM_HC %in% c(0, -1),]
# remove rows with 0 for all human contrasts
orthologs_any_down <- orthologs_any_down[!(orthologs_any_down$sig_MGUS_HC == 0 & orthologs_any_down$sig_SMM_HC == 0 & orthologs_any_down$sig_MM_HC == 0),]

dim(orthologs_any_up) # 18
dim(orthologs_any_down) # 125
# Keep rows only with 0 and 1

orthologs_any_up$gene_name <- gene_anno[match(orthologs_any_up$Human_ensembl_gene_id, gene_anno$gene_id),]$gene_name
orthologs_any_down$gene_name <- gene_anno[match(orthologs_any_down$Human_ensembl_gene_id, gene_anno$gene_id),]$gene_name

sort(unique(orthologs_any_up$gene_name))
sort(unique(orthologs_any_down$gene_name))

Pioneer[(Pioneer %in% orthologs_any_up$gene_name)]
Crisp[(Crisp %in% orthologs_any_up$gene_name)]
Pioneer[(Pioneer %in% orthologs_any_down$gene_name)] #
Crisp[(Crisp %in% orthologs_any_down$gene_name)] #"IRF2"

toPlot_up <- orthologs_any_up
colnames(toPlot_up) <- gsub("sig_", "", colnames(toPlot_up))
library(stringr)
metadata <- as.data.frame(str_split_fixed(colnames(toPlot_up)[3:(ncol(toPlot_up)-1)], "_", n = 2))
colnames(metadata) <- c("Model", "Contrast")
metadata$Model <- gsub("MGUS", "Human", metadata$Model)
metadata$Model <- gsub("SMM", "Human", metadata$Model)
metadata$Model <- gsub("MM", "Human", metadata$Model)
metadata$Contrast[9:11] <- c("MGUS", "SMM", "MM")

rownames(metadata) <- colnames(toPlot_up)[3:(ncol(toPlot_up)-1)]


color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")

color_models <- c("#368716", "#CAB2DC",  "#FB9A99", "#CAF2B0", "#1F78B4","yellow2")
names(color_models) <-c ("Control","CyclinD1", "Mmset", "Trp53", "MIc", "Human")


metadata$Model <- factor(metadata$Model, levels = c("CyclinD1", "Mmset", "Trp53", "MIc", "Human"))
metadata$Combined <- paste(metadata$Model, metadata$Contrast, sep = "_")
metadata$Combined <- gsub("Human_MGUS", "MGUS_HC", metadata$Combined)
metadata$Combined <- gsub("Human_SMM", "SMM_HC", metadata$Combined)
metadata$Combined <- gsub("Human_MM", "MM_HC", metadata$Combined)
metadata$Combined <- factor(metadata$Combined, levels = c( "MGUS_HC", "CyclinD1_MGUS", "Mmset_MGUS", "Trp53_MGUS", "MIc_MGUS","SMM_HC","CyclinD1_MM","Mmset_MM", "Trp53_MM", "MIc_MM", "MM_HC"))

library(ComplexHeatmap)
library(circlize)

ha1 <- HeatmapAnnotation(
  Contrast = metadata$Contrast,
  Model = metadata$Model,
  col = list(Contrast = as.factor(color_stage),
        Model = as.factor(color_models)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


count_table2 <- toPlot_up[,3:(ncol(toPlot_up)-1)]
rownames(count_table2) <- toPlot_up$gene_name
dim(count_table2)
# Order columns to match metadata$Combined
column_ordered <- metadata$Combined[match(colnames(count_table2), rownames(metadata))]

library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs85_any_heatmap.pdf", width=6, height=4)
# order samples by model
Heatmap(as.matrix(count_table2), name = "Differential Expression", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), column_order = order(column_ordered),
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(count_table2)*unit(5, "mm"), #height = nrow(count_table2)*unit(5, "mm"),
        cluster_columns = FALSE, use_raster=TRUE)

dev.off()


toPlot_down <- orthologs_any_down
colnames(toPlot_down) <- gsub("sig_", "", colnames(toPlot_down))

count_table2 <- toPlot_down[,3:(ncol(toPlot_down)-1)]


# Order columns to match metadata$Combined
column_ordered <- metadata$Combined[match(colnames(count_table2), rownames(metadata))]


library(circlize)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
pdf("orthologs_any_down_heatmap.pdf", width=6, height=10)
# order samples by model
Heatmap(as.matrix(count_table2), name = "Differential Expression", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, row_labels = toPlot_down$gene_name,
        show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 4), column_order = order(column_ordered),
        cluster_rows = TRUE, column_split = factor(metadata$Contrast, levels = c("MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"), width = ncol(count_table2)*unit(5, "mm"), #height = nrow(count_table2)*unit(5, "mm"),
        cluster_columns = FALSE, use_raster=TRUE)

dev.off()




# Other genes 
