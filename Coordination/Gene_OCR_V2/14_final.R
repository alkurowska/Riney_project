### Check if annotated genes in differential peaks are differentially expressed 

rm(list = ls())
#################
# LOAD RNA DATA #
#################
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
load("final_genes.RData")
MGUS_HC_genes <- MGUS_HC_filter
SMM_HC_genes <- SMM_HC_filter
MM_HC_genes <- MM_HC_filter

##################
# LOAD ATAC DATA #
##################
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
peak_annotation <- readRDS("peak_annotation.RDS")

# DEA RESULTS - filtered by ChIP
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")
MGUS_HC_atac <- MGUS_HC
SMM_HC_atac <- SMM_HC
MM_HC_atac <- MM_HC

#################
### PROMOTERS ###
#################

# For promoter OCRs, we directly annotated OCRs to genes based on TSS distance (<1kb)

promoters <- peak_annotation[peak_annotation$annotation2 == "Promoter",]
dim(promoters) # 14747 OCRs mapping in promoter regions
# remove rows with NA genes
promoters <- promoters[!is.na(promoters$ENSEMBL),]
# How many OCR-gene pairs ?
dim(promoters) # 13127 OCRs mapping in promoter regions 
summary(is.na(promoters$name))

#################
### ENHANCERS ###
#################

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/GRN/General/enhancers")
enhancers <- read.table("OCR_gene_association.txt", header = TRUE, sep = "\t")
# keep only significant associations p.adj < 0.05
enhancers <- enhancers[enhancers$spearman_adj.pval < 0.05,]
dim(enhancers) # # 39421 OCRs mapping in enhancer regions 
enhancers <- enhancers[!is.na(enhancers$gene_ID),]
# How many OCR-gene pairs ?
dim(enhancers) # # 38549 OCRs mapping in enhancer regions 
summary(is.na(enhancers$peak_ID))

# Are there any peaks in promoters and enhancers
summary(enhancers$peak_ID %in% promoters$name) # NO


# Create OCR_gene pairs for promoters and enhancers
promoters <- promoters[,c("name", "ENSEMBL")]
colnames(promoters) <- c("peak_ID", "gene_ID")
enhancers <- enhancers[,c("peak_ID", "gene_ID")]

# Merge OCR_gene pairs
OCR_gene <- rbind(promoters, enhancers)
dim(OCR_gene) # 51676 OCR-gene pairs

OCR_gene$pairs <- paste0(OCR_gene$peak_ID, "-", OCR_gene$gene_ID)
summary(duplicated(OCR_gene$pairs)) # NO


# Keep only OCR_gene pairs for differential peaks and genes in any of the contrasts 
all_ATAC <- unique(c(MGUS_HC_atac, SMM_HC_atac, MM_HC_atac))
all_genes <- unique(c(MGUS_HC_genes, SMM_HC_genes, MM_HC_genes))

OCR_gene_final <- OCR_gene[OCR_gene$peak_ID %in% all_ATAC & OCR_gene$gene_ID %in% all_genes,]
dim(OCR_gene_final) # 2000 OCR-gene pairs
length(unique(OCR_gene_final$peak_ID)) # 1624 unique peaks
length(unique(OCR_gene_final$gene_ID)) # 1329 unique genes


# save final pairs
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR_V2")
save(OCR_gene_final, file = "ocr_gene_pairs.RData")

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
atac_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
atac_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")


# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
rna_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
rna_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")


OCR_gene_final$MGUS_HC_atac <- 0
OCR_gene_final$SMM_HC_atac <- 0
OCR_gene_final$MM_HC_atac <- 0
OCR_gene_final$MGUS_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_MGUS_HC"]
OCR_gene_final$SMM_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_results)] <- atac_results[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_results)], "sig_SMM_HC"]
OCR_gene_final$MM_HC_atac[OCR_gene_final$peak_ID %in% rownames(atac_res_noscore)] <- atac_res_noscore[OCR_gene_final$peak_ID[OCR_gene_final$peak_ID %in% rownames(atac_res_noscore)], "sig_MM_HC"]

OCR_gene_final$MGUS_HC_rna <- 0
OCR_gene_final$SMM_HC_rna <- 0
OCR_gene_final$MM_HC_rna <- 0
OCR_gene_final$MGUS_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_results)] <- rna_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_results)], "sig_MGUS_HC"]
OCR_gene_final$SMM_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_results)] <- rna_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_results)], "sig_SMM_HC"]
OCR_gene_final$MM_HC_rna[OCR_gene_final$gene_ID %in% rownames(rna_res_noscore)] <- rna_res_noscore[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(rna_res_noscore)], "sig_MM_HC"]

head(OCR_gene_final)



# Check the distribution in contrasts
# DA in MGUS and DE in MGUS 
MGUS_MGUS <- OCR_gene_final[OCR_gene_final$peak_ID %in% MGUS_HC_atac & OCR_gene_final$gene_ID %in% MGUS_HC_genes,]
dim(MGUS_MGUS) # 246 pairs
MGUS_SMM <- OCR_gene_final[OCR_gene_final$peak_ID %in% MGUS_HC_atac & OCR_gene_final$gene_ID %in% SMM_HC_genes,]
dim(MGUS_SMM) # 421 pairs
MGUS_MM <- OCR_gene_final[OCR_gene_final$peak_ID %in% MGUS_HC_atac & OCR_gene_final$gene_ID %in% MM_HC_genes,]
dim(MGUS_MM) # 703 pairs
SMM_SMM <- OCR_gene_final[OCR_gene_final$peak_ID %in% SMM_HC_atac & OCR_gene_final$gene_ID %in% SMM_HC_genes,]
dim(SMM_SMM) # 617 pairs
SMM_MM <- OCR_gene_final[OCR_gene_final$peak_ID %in% SMM_HC_atac & OCR_gene_final$gene_ID %in% MM_HC_genes,]
dim(SMM_MM) # 1098 pairs
MM_MM <- OCR_gene_final[OCR_gene_final$peak_ID %in% MM_HC_atac & OCR_gene_final$gene_ID %in% MM_HC_genes,]
dim(MM_MM) # 1879 pairs

# !!!!! 
MM_SMM <- OCR_gene_final[OCR_gene_final$peak_ID %in% MM_HC_atac & OCR_gene_final$gene_ID %in% SMM_HC_genes,]
dim(MM_SMM) # 937 pairs
MM_MGUS <- OCR_gene_final[OCR_gene_final$peak_ID %in% MM_HC_atac & OCR_gene_final$gene_ID %in% MGUS_HC_genes,]
dim(MM_MGUS) # 468 pairs
SMM_MGUS <- OCR_gene_final[OCR_gene_final$peak_ID %in% SMM_HC_atac & OCR_gene_final$gene_ID %in% MGUS_HC_genes,]
dim(SMM_MGUS) # 326 pairs


length(unique(c(MGUS_MGUS$pairs, SMM_SMM$pairs, MM_MM$pairs))) # 1948

summary(MGUS_MGUS$pairs%in%SMM_SMM$pairs) # 211
summary(MGUS_MGUS$pairs%in%MM_MM$pairs) # 221
summary(SMM_SMM$pairs%in%MM_MM$pairs) # 558

# Is 211 part of 221 ?
MGUS_MGUS_SMM_SMM <- MGUS_MGUS[MGUS_MGUS$pairs%in%SMM_SMM$pairs,]
MGUS_MGUS_MM_MM <- MGUS_MGUS[MGUS_MGUS$pairs%in%MM_MM$pairs,]
summary(MGUS_MGUS_SMM_SMM$pairs%in%MGUS_MGUS_MM_MM$pairs) # 196 TRUE
SMM_SMM_MM_MM <- SMM_SMM[SMM_SMM$pairs%in%MM_MM$pairs,]
summary(MGUS_MGUS_SMM_SMM$pairs%in%SMM_SMM_MM_MM$pairs) # 196 TRUE
MGUS_MGUS_SMM_SMM_MM_MM <- MGUS_MGUS_SMM_SMM[MGUS_MGUS_SMM_SMM$pairs%in%MM_MM$pairs,]

# Keep the extra pairs that are not per contrast
extra_pairs <- OCR_gene_final[!OCR_gene_final$pairs %in% c(MGUS_MGUS$pairs, SMM_SMM$pairs, MM_MM$pairs),]
dim(extra_pairs) # 52 pairs
extra_pairs[,4:9]

# Plot the heatmap 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Create a matrix with contrasts

color_stage <- c("#368716","#D27D46", "#D92F5E", "#8D4E85") 
names(color_stage) <- c("HC", "MGUS_HC", "SMM_HC", "MM_HC")

omic <- c("#166586", "#167F86")
names(omic) <- c("ATAC-Seq", "RNA-Seq")

# create anno dataframe 
anno <- as.data.frame(matrix(NA, nrow = 6, ncol = 2))
rownames(anno) <- colnames(OCR_gene_final)[4:9]
colnames(anno) <- c("Transition", "Omics")

anno$Transition <- c("MGUS_HC", "SMM_HC", "MM_HC", "MGUS_HC", "SMM_HC", "MM_HC")
anno$Omics <- c("ATAC-Seq", "ATAC-Seq", "ATAC-Seq", "RNA-Seq", "RNA-Seq", "RNA-Seq")

ha1 <- rowAnnotation(
  Transition = anno$Transition, 
  col = list(Transition = as.factor(color_stage)),
  na_col = "#D9D9D9", gp = gpar(fontsize = 6),
  annotation_name_gp = gpar(fontsize = 6),   # annotation label size    
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 4),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm")
  ))


count_table2 <- OCR_gene_final[,c("MGUS_HC_atac", "SMM_HC_atac", "MM_HC_atac", "MGUS_HC_rna", "SMM_HC_rna", "MM_HC_rna")]
dim(count_table2) # 2000 x 6

data_to_plot <- count_table2
rownames(data_to_plot) <- OCR_gene_final$pairs

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR_V2")
pdf("atac_rna_heatmap.pdf", width=7, height=3)
ht <- Heatmap(as.matrix(t(data_to_plot)), name = "Differential Signal", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        left_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, row_split = factor(anno$Omics, levels = c("ATAC-Seq", "RNA-Seq")),
        cluster_row_slices = FALSE, column_gap = unit(0, "mm"), row_title_rot = 90,
        cluster_column_slices = FALSE, row_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE, row_title_gp = gpar(fontsize = 6),
        heatmap_legend_param = list( direction = "horizontal",
                title_gp = gpar(fontsize = 6),    # legend title size
                labels_gp = gpar(fontsize = 4),    # legend label size
                grid_height = unit(4, "mm"),       # legend key height
                grid_width = unit(4, "mm")         # legend key width
              ))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# label rows based on ocr-gene contrasts 
all <- intersect(intersect(MGUS_MGUS$pairs, SMM_SMM$pairs), MM_MM$pairs) # 196
length(all)
extra <- extra_pairs$pairs
MGUS_only <- setdiff(MGUS_MGUS$pairs, c(SMM_SMM$pairs, MM_MM$pairs)) # 10
SMM_only <- setdiff(SMM_SMM$pairs, c(MGUS_MGUS$pairs, MM_MM$pairs)) # 44
MM_only <- setdiff(MM_MM$pairs, c(MGUS_MGUS$pairs, SMM_SMM$pairs)) # 1296
MGUS_SMM <- setdiff(intersect(MGUS_MGUS$pairs, SMM_SMM$pairs),all) #15
MGUS_MM <- setdiff(intersect(MGUS_MGUS$pairs, MM_MM$pairs),all) # 25
SMM_MM <- setdiff(intersect(SMM_SMM$pairs, MM_MM$pairs),all) # 362

196+10+44+1296+15+25+362 # 1948

data_to_plot$category <- "MM unique"
data_to_plot$category[rownames(data_to_plot) %in% all] <- "All common"
data_to_plot$category[rownames(data_to_plot) %in% extra] <- "Extra"
data_to_plot$category[rownames(data_to_plot) %in% MGUS_only] <- "MGUS unique"
data_to_plot$category[rownames(data_to_plot) %in% SMM_only] <- "SMM unique"
data_to_plot$category[rownames(data_to_plot) %in% MM_only] <- "MM unique"
data_to_plot$category[rownames(data_to_plot) %in% MGUS_SMM] <- "MGUS_SMM unique"
data_to_plot$category[rownames(data_to_plot) %in% MGUS_MM] <- "MGUS_MM unique"
data_to_plot$category[rownames(data_to_plot) %in% SMM_MM] <- "SMM_MM unique"


summary(as.factor(data_to_plot$category))
pdf("atac_rna_heatmap_labels.pdf", width=10, height=3)
ht <- Heatmap(as.matrix(t(data_to_plot[,1:6])), name = "Differential Signal", col = colorRamp2(c(-1, 0, 1), c("#4475B3", "white", "#D7342A")),
        left_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, 
        # add border of each cell
        #rect_gp = gpar(col = "white", lwd = 0.1),
        # row labels size 
        border = TRUE,
        row_title_gp = gpar(fontsize = 4),
        column_title_gp = gpar(fontsize = 4),
        cluster_rows = TRUE, column_split = factor(data_to_plot$category, levels = c("All common", "MGUS unique", "SMM unique", "MM unique", "MGUS_SMM unique", "MGUS_MM unique", "SMM_MM unique", "Extra")),
        row_split = factor(anno$Omics, levels = c("ATAC-Seq", "RNA-Seq")),
        cluster_row_slices = FALSE, column_gap = unit(0.1, "mm"), row_title_rot = 90,
        cluster_column_slices = FALSE, row_gap = unit(0.5, "mm"), column_title_rot = 90,
        cluster_columns = TRUE, use_raster=TRUE, column_title_side = "bottom",
        heatmap_legend_param = list(
                title_gp = gpar(fontsize = 6),    # legend title size
                labels_gp = gpar(fontsize = 4),    # legend label size
                grid_height = unit(4, "mm"),       # legend key height
                grid_width = unit(4, "mm")         # legend key width
              ))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# summarize the number of genes in each category
dynamics <- data_to_plot[,1:6]
dynamics <- dynamics[!rownames(dynamics) %in% extra,]
dim(dynamics) # 1441 x 6
library(stringr)

MGUS_up <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$MGUS_HC_rna == 1,]),"-",2)) # 114
SMM_up <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$SMM_HC_atac == 1 & dynamics$SMM_HC_rna == 1,]),"-",2)) # 174
MM_up <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MM_HC_atac == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 387

MGUS_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$MGUS_HC_rna == -1,]),"-",2)) # 78
SMM_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$SMM_HC_atac == -1 & dynamics$SMM_HC_rna == -1,]),"-",2)) # 335
MM_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MM_HC_atac == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 1033

colnames(MGUS_up) <- c("peak_ID", "gene_ID")
colnames(SMM_up) <- c("peak_ID", "gene_ID")
colnames(MM_up) <- c("peak_ID", "gene_ID")
colnames(MGUS_down) <- c("peak_ID", "gene_ID")
colnames(SMM_down) <- c("peak_ID", "gene_ID")
colnames(MM_down) <- c("peak_ID", "gene_ID")

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
saveRDS(MGUS_up, file = "MGUS_up.rds")
saveRDS(SMM_up, file = "SMM_up.rds")
saveRDS(MM_up, file = "MM_up.rds")
saveRDS(MGUS_down, file = "MGUS_down.rds")
saveRDS(SMM_down, file = "SMM_down.rds")
saveRDS(MM_down, file = "MM_down.rds")


MGUS_atac_up_rna_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$MGUS_HC_rna == -1,]),"-",2)) # 14
SMM_atac_up_rna_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$SMM_HC_atac == 1 & dynamics$SMM_HC_rna == -1,]),"-",2)) # 16
MM_atac_up_rna_down <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MM_HC_atac == 1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 30

colnames(MGUS_atac_up_rna_down) <- c("peak_ID", "gene_ID")
colnames(SMM_atac_up_rna_down) <- c("peak_ID", "gene_ID")
colnames(MM_atac_up_rna_down) <- c("peak_ID", "gene_ID")

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
saveRDS(MGUS_atac_up_rna_down, file = "MGUS_atac_up_rna_down.rds")
saveRDS(SMM_atac_up_rna_down, file = "SMM_atac_up_rna_down.rds")
saveRDS(MM_atac_up_rna_down, file = "MM_atac_up_rna_down.rds")

# ATAC MGUS open, SMM open, MM open -> RNA MGUS up, SMM up, MM up
group1 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 1 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 88

# ATAC MGUS open, SMM open, MM open -> RNA MGUS non, SMM up, MM up
group2 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 19

# ATAC MGUS open, SMM open, MM open -> RNA MGUS non, SMM non, MM up
group3 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 44

# ATAC MGUS non, SMM open, MM open -> RNA MGUS non, SMM up, MM up
group4 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 16

# ATAC MGUS non, SMM open, MM open -> RNA MGUS non, SMM non, MM up
group5 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 37

# ATAC MGUS non, SMM non, MM open -> RNA MGUS non, SMM non, MM up
group6 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 0 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 73



# Same for down
# ATAC MGUS down, SMM down, MM down -> RNA MGUS down, SMM down, MM down
group7 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == -1 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 78

# ATAC MGUS down, SMM down, MM down -> RNA MGUS non, SMM down, MM down,"-",2)
group8 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 119

# ATAC MGUS down, SMM down, MM down -> RNA MGUS non, SMM non, MM down
group9 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 138

# ATAC MGUS non, SMM down, MM down -> RNA MGUS non, SMM down, MM down
group10 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 77

# ATAC MGUS non, SMM down, MM down -> RNA MGUS non, SMM non, MM down
group11 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 126

# ATAC MGUS non, SMM non, MM down -> RNA MGUS non, SMM non, MM down
group12 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 0 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 256



# Mixed ATAC closing -> RNA up
# ATAC MGUS down, SMM down, MM down -> RNA MGUS up, SMM up, MM up
group19 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 1 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 16

# ATAC MGUS down, SMM down, MM down -> RNA MGUS non, SMM up, MM up
group20 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 14

# ATAC MGUS down, SMM down, MM down -> RNA MGUS non, SMM non, MM up
group21 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == -1 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 57

# ATAC MGUS non, SMM down, MM down -> RNA MGUS non, SMM up, MM up
group22 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 1 & dynamics$MM_HC_rna == 1,]),"-",2)) # 14

# ATAC MGUS non, SMM down, MM down -> RNA MGUS non, SMM non, MM up
group23 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == -1 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 43

# ATAC MGUS non, SMM non, MM down -> RNA MGUS non, SMM non, MM up
group24 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 0 & dynamics$MM_HC_atac == -1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == 1,]),"-",2)) # 95




# ATAC opening -> RNA down

# ATAC MGUS open, SMM open, MM open -> RNA MGUS down, SMM down, MM down
group13 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == -1 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 14

# ATAC MGUS open, SMM open, MM open -> RNA MGUS non, SMM down, MM down
group14 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 16

# ATAC MGUS open, SMM open, MM open -> RNA MGUS non, SMM non, MM down
group15 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 1 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 30

# ATAC MGUS non, SMM open, MM open -> RNA MGUS non, SMM down, MM down
group16 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == -1 & dynamics$MM_HC_rna == -1,]),"-",2)) # 7

# ATAC MGUS non, SMM open, MM open -> RNA MGUS non, SMM non, MM down
group17 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 1 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 28

# ATAC MGUS non, SMM non, MM open -> RNA MGUS non, SMM non, MM down
group18 <- as.data.frame(str_split_fixed(rownames(dynamics[dynamics$MGUS_HC_atac == 0 & dynamics$SMM_HC_atac == 0 & dynamics$MM_HC_atac == 1 & dynamics$MGUS_HC_rna == 0 & dynamics$SMM_HC_rna == 0 & dynamics$MM_HC_rna == -1,]),"-",2)) # 36


# sum all
sum(c(88, 19, 44, 16, 37, 73, 78, 119, 138, 77, 126, 256, 16, 14, 57, 14, 43, 95, 14, 16, 30, 7, 28, 36)) # 1441

# 1948 - 1441 = 507


# PLOT sankey

#### PREP THE DATA ####
colnames(dynamics)
toPlot <- dynamics[,c("MGUS_HC_atac", "MGUS_HC_rna", "SMM_HC_atac", "SMM_HC_rna", "MM_HC_atac", "MM_HC_rna")]

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Pair <- rownames(dynamics)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial

genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(MGUS_HC_atac, MGUS_HC_rna, SMM_HC_atac, SMM_HC_rna, MM_HC_atac, MM_HC_rna),
  names_to = "Condition",
  values_to = "Category"
)

genes_long$Group <- gsub("_(atac|rna)", "", genes_long$Condition)
genes_long$Modality <- gsub(".*_(atac|rna)", "\\1", genes_long$Condition)

genes_long$Group <- factor(genes_long$Group, levels = c("MGUS_HC", "SMM_HC", "MM_HC"))
genes_long$Modality <- factor(genes_long$Modality, levels = c("atac", "rna"))

genes_long <- genes_long[order(genes_long$Group, genes_long$Modality), ]
genes_long$Condition <- interaction(genes_long$Group, genes_long$Modality, sep = "_")

genes_long$Condition <- factor(genes_long$Condition, levels = unique(genes_long$Condition))

# Add spacer levels
genes_long$Display_Condition <- as.character(genes_long$Condition)
genes_long$Display_Condition <- factor(genes_long$Display_Condition,
  levels = c(
    "MGUS_HC_atac", "MGUS_HC_rna", "spacer1",
    "SMM_HC_atac", "SMM_HC_rna", "spacer2",
    "MM_HC_atac", "MM_HC_rna"
  )
)

# Plot
getwd()
png("ocr_gene_sanky.png", width=2000, height=1000, res=300)

# Plot with spacer strata dropped
ggplot(genes_long[!grepl("spacer", genes_long$Display_Condition), ], aes(
  x = Display_Condition, stratum = Category, alluvium = Pair,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium", alpha = 0.8, aes(order = as.integer(Display_Condition))) +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "OCR-Gene Pairs Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", 
                               "Down-regulated" = "#4475B3", 
                               "Non-differential" = "#F0F0F0")) +
  scale_x_discrete(drop = FALSE)
dev.off()





colnames(group1) <- c("peak_ID", "gene_ID")
colnames(group2) <- c("peak_ID", "gene_ID")
colnames(group3) <- c("peak_ID", "gene_ID")
colnames(group4) <- c("peak_ID", "gene_ID")
colnames(group5) <- c("peak_ID", "gene_ID")
colnames(group6) <- c("peak_ID", "gene_ID")
colnames(group7) <- c("peak_ID", "gene_ID")
colnames(group8) <- c("peak_ID", "gene_ID")
colnames(group9) <- c("peak_ID", "gene_ID")
colnames(group10) <- c("peak_ID", "gene_ID")
colnames(group11) <- c("peak_ID", "gene_ID")
colnames(group12) <- c("peak_ID", "gene_ID")
colnames(group13) <- c("peak_ID", "gene_ID")
colnames(group14) <- c("peak_ID", "gene_ID")
colnames(group15) <- c("peak_ID", "gene_ID")
colnames(group16) <- c("peak_ID", "gene_ID")
colnames(group17) <- c("peak_ID", "gene_ID")
colnames(group18) <- c("peak_ID", "gene_ID")
colnames(group19) <- c("peak_ID", "gene_ID")
colnames(group20) <- c("peak_ID", "gene_ID")
colnames(group21) <- c("peak_ID", "gene_ID")
colnames(group22) <- c("peak_ID", "gene_ID")
colnames(group23) <- c("peak_ID", "gene_ID")
colnames(group24) <- c("peak_ID", "gene_ID")


setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
save(group1, group2, group3, group4, group5, group6, 
     group7, group8, group9, group10, group11, group12,
     group13, group14, group15, group16, group17, group18,
     group19, group20, group21, group22, group23, group24, file = "OCR-gene_pairs_GROUPS.RData")




