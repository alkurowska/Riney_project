#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: chromVAR_TF.R          #
#***********************************#


# https://greenleaflab.github.io/chromVAR/articles/Introduction.html
# https://bioconductor.org/packages/3.3/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment

# loading required packages
# in bash: mamba activate chromVAR
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(GenomicRanges)
library(ggplot2)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(purrr)
library(TFBSTools)
library(JASPAR2024)
library(ComplexHeatmap)


## Getting the required files
# 1- Set of OCR
# Consensus Peaks for the final significant peaks
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_peaks.RData")

final_peaks <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# 2- OCR counts and metadata
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata")
my_counts_matrix <- Counts$counts

peaks <- Counts$annotation
rownames(peaks) <- paste(peaks$Chr, peaks$Start, peaks$End, sep = "_")
rownames(my_counts_matrix) <- rownames(peaks)

# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]

# change colnames of the count data
# Remove everything after "_S"
colnames(my_counts_matrix) <- gsub("_S.*", "", colnames(my_counts_matrix))

# Load metadata

# load metadata
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
fish <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(fish) <- fish[1,]
fish <- fish[-1,]

# Remove FISH data with NA values
fish <- fish[!is.na(fish$`1q amp`),] # 202
dim(fish)

# Remove FISH data with empty RNA seq sampels
fish <- fish[!fish$`Sample RNAseq` == "",] # 197 x 17

my_counts_matrix <- my_counts_matrix[,colnames(my_counts_matrix)%in%fish$`Sample ATACseq`]
dim(my_counts_matrix) # 157817  x  185

fish <- fish[fish$`Sample ATACseq` %in% colnames(my_counts_matrix),]
dim(fish) # 185 x 17

rownames(fish) <- fish$`Sample ATACseq`

# Samples/peaks to keep
my_counts_matrix <- my_counts_matrix[,rownames(fish)]
dim(my_counts_matrix) # 157817  x  185

### Filtered peaks
### STEP 1: LOAD DATA
# Create the GRanges of consensus peaks
peaks <- peaks[final_peaks,]
peak_granges <- makeGRangesFromDataFrame(peaks, keep.extra.columns = T)
peak_granges

my_counts_matrix <- my_counts_matrix[final_peaks,]
table(rownames(my_counts_matrix) == names(peak_granges))



## Creating the object 
table(rownames(fish) == colnames(my_counts_matrix))
fragment_counts <- SummarizedExperiment(assays=list(counts=my_counts_matrix),
                                        rowRanges=peak_granges, colData=fish)

# Getting GC content of peaks

fragment_counts2 <- addGCBias(fragment_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg38)

#Filtering peaks
counts_filtered <- filterPeaks(fragment_counts2, non_overlapping = TRUE)


#Get motifs and what peaks contain motifs
## JASPAR 2024
#https://academic.oup.com/nar/article/46/D1/D252/4616875

#Read the motifs
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
motifs24 <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))

# Save motifs names 

# Save motifs names 
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
load("JASPAR2024.RData")


# keep only expressed
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
load("norm_to_voom.RData")
genes <- unique(gene_anno[gene_anno$gene_id%in%rownames(d1_norm$counts),]$gene_name)
TFs <- unique(names[,1][names[,1]%in%genes]) #453

motifs24 <- motifs24[which(names[,1]%in%TFs)] #484

# Get motif matches in peaks
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/03_chromVAR")
motif_ix <- matchMotifs(motifs24, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#Save as RData
save(motif_ix, file = "Jaspar_motifs_matches.RData")                        
# Get motif positions within peaks for motifs in peaks 

jaspar_ranges <- matchMotifs(motifs24, peak_granges, genome = "hg38", out = "positions") 
# for(i in 1:length(motifs24)){
#   names(jaspar_ranges)[i] <- unname(name(motifs24[i]))
# }
save(jaspar_ranges, file = "Jaspar_motifs_positions.RData")
##>>>>>>>>>>>>>>>>>>>>> Relate each OCR to motif

jaspar_names <- names(jaspar_ranges)

ocr_jaspar_matrix <- matrix(0,nrow=nrow(peaks),ncol=length(jaspar_names))
colnames(ocr_jaspar_matrix) <- jaspar_names
names(peak_granges) <- paste0(peaks$Chr,"_", peaks$Start, "_", peaks$End)
rownames(ocr_jaspar_matrix) <- names(peak_granges)

for(i in 1:length(jaspar_names)){
  cat(i," - ")
  #select granges for target motif
  jaspar_gr_motif <- jaspar_ranges[[i]]
  #ranges overlaping
  ocr_jaspar_association <- findOverlaps(jaspar_gr_motif,peak_granges,type = "within")
  if(length(ocr_jaspar_association)>0){
    ocr_jaspar_matrix[unique(ocr_jaspar_association@to),i] <- 1    
  }  
}

colnames(ocr_jaspar_matrix) <- jaspar_names

summary(colSums(ocr_jaspar_matrix)==0) #FALSE
write.table(ocr_jaspar_matrix, "ocr_jaspar_matrix.txt", sep = "\t", col.names = T, row.names = T)


summary(rowSums(ocr_jaspar_matrix)==0) #FALSE

#Get annotations
anno_jaspar <- getAnnotations(ocr_jaspar_matrix[rownames(counts_filtered),], 
                         rowRanges = rowRanges(counts_filtered))


######### computing deviations
#Background Peaks
bg <- getBackgroundPeaks(object = counts_filtered)

dev_jaspar <- computeDeviations(object = counts_filtered, background = bg,
                         annotations = motif_ix)

save(list=c("dev_jaspar"),file=paste(format(Sys.Date(),"%Y%m%d"),"_","dev_jaspar","_dev.RData",sep=""))

# variability
variability <- computeVariability(dev_jaspar)
toKeep <- names[rownames(variability),]


cat("Top variable JASPAR motifs:\n")
print(variability[order(variability$variability, decreasing = TRUE)[1:20],])
png(paste(format(Sys.Date(),"%Y%m%d"),"_","JASPAR_variability.png",sep=""))
plotVariability(variability, use_plotly = FALSE)
dev.off()

save(list=c("variability"),file=paste(format(Sys.Date(),"%Y%m%d"),"_","JASPAR","_variability.RData",sep=""))
write.table(variability, "variability_jaspar.txt", sep = "\t")


png("histogram_variability.png")
hist(variability$variability)
dev.off()
quantile(variability$variability,0.25) # cut-off: 0.9027363

# clustering samples


### STEP 2: Multiple model analysis
fish_atac <- fish

colnames(fish_atac)[4] <- "del17p"
colnames(fish_atac)[5] <- "del1p"
colnames(fish_atac)[6] <- "trans_4_14"
colnames(fish_atac)[7] <- "trans_14_16"
colnames(fish_atac)[11] <- "qall"
fish_atac$qall <- factor(fish_atac$qall)
fish_atac$del1p <- factor(fish_atac$del1p)
fish_atac$del17p <- factor(fish_atac$del17p)
fish_atac$trans_4_14 <- factor(fish_atac$trans_4_14)
fish_atac$trans_14_16 <- factor(fish_atac$trans_14_16)

fish_anno <- fish_atac
fish_anno$translocation <- as.character(fish_anno$trans_4_14)
fish_anno[fish_anno$trans_14_16 == "t(14;16)",]$translocation <- "t(14;16)"
fish_anno$translocation <- as.character(fish_anno$translocation)

fish_anno$qall <- as.character(fish_anno$qall)
fish_anno$del1p <- as.character(fish_anno$del1p)
fish_anno[fish_anno$del1p == "1p del",]$del1p <- "del(1p)"
fish_anno$del17p <- as.character(fish_anno$del17p)
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

#plot

library(circlize)
f1 <- colorRamp2(c(-1, 1), c("blue", "red"))

sample_cor_ja <- getSampleCorrelation(dev_jaspar)

png(paste(format(Sys.Date(),"%Y%m%d"),"_","JASPAR_sample_correlation.png",sep=""), width=5000, height=5000, res=300)
draw(Heatmap(sample_cor_ja, top_annotation = ha1, 
             show_row_names = FALSE, show_column_names =FALSE, row_labels = NULL,  cluster_column_slices = FALSE,
              cluster_row_slices = FALSE,
             row_names_gp = gpar(fontsize = 3), col=f1, column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
             name="rho", row_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
             column_title = "Correlation for bias corrected deviations"))
dev.off()


# Get deviations and scores
#bias corrected deviations in accessibility

TF_deviations_ja <- deviations(dev_jaspar) 
write.table(TF_deviations_ja, "TF_deviations_jaspar.txt", sep = "\t")


#deviationScores are the Z-scores for each bias corrected deviations
TF_scores_ja <- deviationScores(dev_jaspar)
write.table(TF_scores_ja, "jaspar_TFBS_score.txt", sep = "\t")


# Graphical exploration
#heatmap
#JASPAR
variability_jaspar <- variability
target_TF_ja <- rownames(variability_jaspar[order(variability_jaspar$variability, decreasing = TRUE)[1:20],])

target_TF_position_ja <- vector()

for(i in 1:length(target_TF_ja)){
  target_TF_position_ja <- c(target_TF_position_ja,grep(target_TF_ja[i],rownames(TF_scores_ja)))
}

toLabel <- rownames(TF_scores_ja)[target_TF_position_ja]
toLabel <- name(motifs24[toLabel])

rha = rowAnnotation(foo = anno_mark(at = target_TF_position_ja, 
                                    labels = toLabel,
                                    labels_gp = gpar(fontsize = 6)))

f1 <- colorRamp2(c(-1, 1), c("darkviolet", "yellow"))
##>>>>>>>>>>>>>>>>>>>>> Heatmap

jpeg("TFBS_accessibility_JASPAR.jpg", width=5000, height=5000, res=600)
Heatmap(TF_scores_ja, top_annotation = ha1, right_annotation = rha, 
        show_row_names = FALSE, show_column_names =FALSE, col = f1,
        column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        column_title_gp = gpar(fontsize = 10),
        cluster_columns = T, 
        cluster_column_slices = FALSE,
        cluster_rows = T,
        row_names_gp = gpar(fontsize = 3), 
        name="Z-score",
        column_title = "TFBS accessibility score")

dev.off()


f1 <- colorRamp2(c(-0.1, 0, 0.1), c("purple", "white" , "yellow"))
jpeg("Bias-corrected_deviations_JASPAR.jpg", width=5000, height=5000, res=600)
Heatmap(as.matrix(TF_deviations_ja), top_annotation = ha1, right_annotation = rha, 
        show_row_names = FALSE, show_column_names =FALSE, col = f1,
        column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        column_title_gp = gpar(fontsize = 10),
        cluster_columns = T, 
        cluster_column_slices = FALSE,
        cluster_rows = T,
        row_names_gp = gpar(fontsize = 3), 
        name="Bias
correct
devs",
        column_title = "TFBS accessibility deviations")

dev.off()




# # Differential analysis
# # The differentialDeviations function determines whether there is a significant difference between the bias corrected deviations for a given annotation between different groups. The groups can be specified by giving the column name of a column in the colData of the dev object or a vector of group assignments.

# diff_acc <- differentialDeviations(dev_jaspar, "Stage")


# #The differentialVariability function determines whether there is a significant difference between the variability of any of the annotations between different groups.

# diff_var <- differentialVariability(dev_jaspar, "Stage")







# Limma analysis
load("20250328_dev_jaspar_dev.RData")

# First, the normalization of z-scores
assays(dev_jaspar)$norm <- scale(assays(dev_jaspar)$z)
#devMat <- assays(dev_jaspar)$norm

devMat <- assays(dev_jaspar)$z
  

# Estimate the fold changes and standard errors by fitting a linear model for each gene
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("norm_to_voom.RData")

#design
fish_atac <- fish_atac
colnames(fish_atac)[4] <- "del17p"
colnames(fish_atac)[5] <- "del1p"
colnames(fish_atac)[6] <- "trans_4_14"
colnames(fish_atac)[7] <- "trans_14_16"
colnames(fish_atac)[11] <- "qall"

fish_atac$qall <- gsub("neutral", "", fish_atac$qall)
fish_atac$del17p <- gsub("neutral", "", fish_atac$del17p)
fish_atac$del1p <- gsub("neutral", "", fish_atac$del1p)
fish_atac$trans_4_14 <- gsub("neutral", "", fish_atac$trans_4_14)
fish_atac$trans_14_16 <- gsub("neutral", "", fish_atac$trans_14_16)

# Relevant variables for the design matrix
stage = as.factor(fish_atac$Stage)
amp = as.factor(fish_atac$qall)
del1p = as.factor(fish_atac$del1p)
del17p = as.factor(fish_atac$del17p)
trans_14_16 = as.factor(fish_atac$trans_14_16)
trans_4_14 = as.factor(fish_atac$trans_4_14)
sex = as.factor(fish_atac$Sex)
gsva <- as.numeric(fish_atac$gsva)
frip <- as.numeric(fish_atac$FRIP)

# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16 + gsva + frip)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"

fit <- lmFit(devMat, design)
  
# Apply empirical Bayes smoothing to the standard errors
  
fit2 <- eBayes(fit)
  
# Show statistics for the top 10 genes
summary(decideTests(fit2))
#       Intercept MGUS  MM SMM male qaber del1p del17 trans_4_14 trans_14_16
# Down          14  137 161 140    0     0     0     0          0           0
# NotSig       289  334 304 328  484   484   484   484        484         484
# Up           181   13  19  16    0     0     0     0          0           0
#       gsva frip
#Down      0    0
#NotSig  484  484
#Up        0    0
  
# # Check the difference between samples with 1q amp, 1q gain and 1q gainORamp for SMM and MM samples 
myContrMatrix = makeContrasts(
   MM_SMM = MM - SMM, #trasitiion
   SMM_MGUS = SMM - MGUS, #trasitiion
   MGUS_HC = MGUS , #trasitiion
   MM_MGUS = MM - MGUS, #transition
   MM_HC = MM , #trasitiion
   SMM_HC = SMM , #trasitiion
   levels= colnames(design)
  )


## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(devMat), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #119
sum(contrast_stats$adj.P.Val < 0.01) #77

dtest <- decideTests(fit.contrasts2)

contr.matrix <- myContrMatrix
all_pval <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(devMat))

for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i, number=nrow(devMat), sort.by="none")
  
  all_pval[,i] <- target_stats[,"P.Value"]
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] > 0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] < -0.58))
  sig[sel2,i] <- -1 
}

colnames(all_pval) = paste("p.value_",colnames(contr.matrix),sep="")
colnames(all_padj) = paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) = paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) = paste("sig_",colnames(contr.matrix),sep="")

dea_results = data.frame(sig,all_padj,all_fch,all_pval)
rownames(dea_results) = rownames(contrast_stats)
which(duplicated(rownames(dea_results)))

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0]) #205
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #205
#with a adj.pval<0.05 & |FC|>2 (logFC>0.58)

table(rownames(dea_results) == names(motifs24)) #TRUE
toName <- names[names(motifs24),]
dea_results$TFs <- toName
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/03_chromVAR")
write.table(dea_results,"chromVAR_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential for each contrast
length(which(dea_results$sig_SMM_HC != 0)) # 156
length(which(dea_results$sig_MGUS_HC != 0)) # 150

table(dea_results$sig_SMM_HC)
table(dea_results$sig_MGUS_HC)

length(intersect(which(dea_results$sig_SMM_HC == 1), which(dea_results$sig_MGUS_HC == 1))) # 13
length(intersect(which(dea_results$sig_SMM_HC == -1), which(dea_results$sig_MGUS_HC == -1))) # 111



#### NO SCORE


# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16 + frip)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"

fit <- lmFit(devMat, design)
  
# Apply empirical Bayes smoothing to the standard errors
  
fit2 <- eBayes(fit)
  
# Show statistics for the top 10 genes
summary(decideTests(fit2))
#       Intercept MGUS  MM SMM male qaber del1p del17 trans_4_14 trans_14_16
# Down          25  120 146 126    0     0     0     0          0           0
# NotSig       284  351 305 330  484   484   484   484        484         484
# Up           175   13  33  28    0     0     0     0          0           0
#        frip
# Down      0
# NotSig  484
# Up        0
  
# # Check the difference between samples with 1q amp, 1q gain and 1q gainORamp for SMM and MM samples 
myContrMatrix = makeContrasts(
   MM_SMM = MM - SMM, #trasitiion
   SMM_MGUS = SMM - MGUS, #trasitiion
   MGUS_HC = MGUS , #trasitiion
   MM_MGUS = MM - MGUS, #transition
   MM_HC = MM , #trasitiion
   SMM_HC = SMM , #trasitiion
   levels= colnames(design)
  )


## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(devMat), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #114
sum(contrast_stats$adj.P.Val < 0.01) #83

dtest <- decideTests(fit.contrasts2)

contr.matrix <- myContrMatrix
all_pval <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(devMat))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(devMat))

for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i, number=nrow(devMat), sort.by="none")
  
  all_pval[,i] <- target_stats[,"P.Value"]
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] > 0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] < -0.58))
  sig[sel2,i] <- -1 
}

colnames(all_pval) = paste("p.value_",colnames(contr.matrix),sep="")
colnames(all_padj) = paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) = paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) = paste("sig_",colnames(contr.matrix),sep="")

dea_results = data.frame(sig,all_padj,all_fch,all_pval)
rownames(dea_results) = rownames(contrast_stats)
which(duplicated(rownames(dea_results)))

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0]) #205
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #205
#with a adj.pval<0.05 & |FC|>2 (logFC>0.58)

table(rownames(dea_results) == names(motifs24)) #TRUE
toName <- names[names(motifs24),]
dea_results$TFs <- toName
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/03_chromVAR")
write.table(dea_results,"chromVAR_dea_results_noscore.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential for each contrast
length(which(dea_results$sig_MM_HC != 0)) # 179

table(dea_results$sig_MM_HC)