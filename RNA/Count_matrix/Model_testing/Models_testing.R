###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           Models_testing.R        ####
####              HUMAN                ####
###########################################


# Setting working directory
set.seed(1234)

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("voom_to_models_testing.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
norm_data <- norm_data

# Load Cytogenetics 
fish_rna <- fish_rna
colnames(fish_rna)[4] <- "del17p"
colnames(fish_rna)[5] <- "del1p"
colnames(fish_rna)[6] <- "trans_4_14"
colnames(fish_rna)[7] <- "trans_14_16"
colnames(fish_rna)[11] <- "qall"

fish_rna$del17p <- gsub("neutral", "", fish_rna$del17p)
fish_rna$del1p <- gsub("neutral", "", fish_rna$del1p)
fish_rna$trans_4_14 <- gsub("neutral", "", fish_rna$trans_4_14)
fish_rna$trans_14_16 <- gsub("neutral", "", fish_rna$trans_14_16)
fish_rna$qall <- gsub("neutral", "", fish_rna$qall)

# Relevant covariates 
Stage = as.factor(fish_rna$Stage)
Sex = as.factor(fish_rna$Sex)
qall = as.factor(fish_rna$qall)
del1p = as.factor(fish_rna$del1p)
del17p = as.factor(fish_rna$del17p)
trans_14_16 = as.factor(fish_rna$trans_14_16)
trans_4_14 = as.factor(fish_rna$trans_4_14)


# Add GSVA scores
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
table(rownames(fish_rna) == rownames(scores))
fish_rna$gsva <- scores[rownames(fish_rna), "GSVA_low"]
gsva <- as.numeric(fish_rna$gsva)

library("variancePartition")

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/fish_rna about each sample

# Stage together + 1q together
form <- ~ (1 | Stage) + (1 | Sex) + (1 | qall) + (1 | del1p) + (1 | del17p) + (1 | trans_14_16) + (1 | trans_4_14) + gsva

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, fish_rna)
colnames(C) <- c("Stage", "Sex", "1q aberration", "1p del", "17p del", "t(14;16)", "t(4;14)", "GSVA")
rownames(C) <- colnames(C)  

# Plot correlation matrix
# between all pairs of variables

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
png("correlation_matrix.png", res = 300, width = 2000, height = 2000)
plotCorrMatrix(C)
dev.off()

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,

varPart <- fitExtractVarPartModel(norm_data, form, fish_rna)

#Rturns an object that stores the variance fractions for each gene and each variable in the formula specified. These fractions can be accessed just like a data.frame:


# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
colnames(vp) <- c("GSVA", "Stage","t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")

# Figure 1b
# violin plot of contribution of each variable to total variance
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
png("variance_explained.png", res = 300, width = 4000, height = 2000)
plotVarPart(vp)
dev.off()

#Genes important for aberrations
amp <- gene_anno[gene_anno$gene_name %in% c("PYGO2", "CKS1B", "ZBTB7B"),]$gene_id
del1 <- gene_anno[gene_anno$gene_name %in% c( "FAF1", "CDKN2C"),]$gene_id
del17 <- gene_anno[gene_anno$gene_name %in% c("TP53", "SAT2", "EFNB3" ),]$gene_id
trans <- gene_anno[gene_anno$gene_name %in% c("MAF", "FGFR3"),]$gene_id

# Figure 1a

# Bar plot of variance fractions for the first 10 genes
png("variance_fractions.png", res = 300, width = 6, height = 4, units = "in")
ids <- c(amp, del1, del17, trans)
toPlot <- vp[rownames(vp)%in%ids, ]
rownames(toPlot) <- gene_anno[gene_anno$gene_id %in% rownames(toPlot),]$gene_name
plotPercentBars(toPlot)
dev.off()

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
save(varPart, file = "varPart_data.RData")


# Plot Top 15 genes explained by GSVA
# Sort by GSVA
vp_gsva <- vp[order(vp[,"GSVA"], decreasing = TRUE),]
vp_gsva <- vp_gsva[1:30,]
rownames(vp_gsva) <- gene_anno[gene_anno$gene_id %in% rownames(vp_gsva),]$gene_name
toPlot <- vp_gsva 
png("var_frac_top_gsva.png", res = 300, width = 6, height = 4, units = "in")
plotPercentBars(toPlot)
dev.off()


# check all of the APOBEC genes
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
load("varPart_data.RData")
# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
colnames(vp) <- c("GSVA", "Stage","t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")
toKeep <- gene_anno[grep("APOBEC",gene_anno$gene_name),]$gene_id

vp_apobec <- vp[rownames(vp) %in% toKeep,]
rownames(vp_apobec) <- gene_anno[gene_anno$gene_id %in% rownames(vp_apobec),]$gene_name
toPlot <- vp_apobec
png("var_frac_apobec.png", res = 300, width = 6, height = 4, units = "in")
plotPercentBars(toPlot)
dev.off()



