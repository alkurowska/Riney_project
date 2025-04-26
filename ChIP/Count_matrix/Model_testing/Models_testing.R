###########################################
####     ATAC-Seq pipeline - DOUBLE    ####
####           Models_testing.R        ####
####              HUMAN                ####
###########################################


# Setting working directory
set.seed(1234)

setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("voom_to_models_testing.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
norm_data <- norm_data

# Load Cytogenetics 
fish_chip <- fish_chip
colnames(fish_chip)[4] <- "del17p"
colnames(fish_chip)[5] <- "del1p"
colnames(fish_chip)[6] <- "trans_4_14"
colnames(fish_chip)[7] <- "trans_14_16"
colnames(fish_chip)[11] <- "qall"

fish_chip$del17p <- gsub("neutral", "", fish_chip$del17p)
fish_chip$del1p <- gsub("neutral", "", fish_chip$del1p)
fish_chip$trans_4_14 <- gsub("neutral", "", fish_chip$trans_4_14)
fish_chip$trans_14_16 <- gsub("neutral", "", fish_chip$trans_14_16)
fish_chip$qall <- gsub("neutral", "", fish_chip$qall)

# Relevant covariates 
Stage = as.factor(fish_chip$Stage)
Sex = as.factor(fish_chip$Sex)
qall = as.factor(fish_chip$qall)
del1p = as.factor(fish_chip$del1p)
del17p = as.factor(fish_chip$del17p)
trans_14_16 = as.factor(fish_chip$trans_14_16)
trans_4_14 = as.factor(fish_chip$trans_4_14)


# Add GSVA scores
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
fish_chip <- fish_chip[fish_chip$`Sample RNAseq`%in%rownames(scores),]
dim(fish_chip) #128 x 18
fish_chip$gsva <- scores[fish_chip$`Sample RNAseq`, "GSVA_low"]

library("variancePartition")

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/fish_atac about each sample

# Stage
form <- ~ (1 | Stage) + (1 | Sex) + (1 | qall) + (1 | del1p) + (1 | del17p) + (1 | trans_14_16) + (1 | trans_4_14) + gsva

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, fish_chip)
colnames(C) <- c("Stage", "Sex", "1q aberration", "1p del", "17p del", "t(14;16)", "t(4;14)", "GSVA")
rownames(C) <- colnames(C)  

# Plot correlation matrix
# between all pairs of variables

setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
png("correlation_matrix.png", res = 300, width = 2000, height = 2000)
plotCorrMatrix(C)
dev.off()

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,

varPart <- fitExtractVarPartModel(norm_data, form, fish_chip)

#Rturns an object that stores the variance fractions for each gene and each variable in the formula specified. These fractions can be accessed just like a data.frame:


# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
colnames(vp) <- c("GSVA", "17p del",  "1p del", "1q aber", "Sex", "Stage", "t(14;16)", "t(4;14)", "Residuals")
vp <- vp[,c("GSVA", "Stage", "t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")]

# Figure 1b
# violin plot of contribution of each variable to total variance
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
png("variance_explained.png", res = 300, width = 4000, height = 2000)
plotVarPart(vp)
dev.off()

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
save(varPart, file = "varPart_data.RData")



# No score
form <- ~ (1 | Stage) + (1 | Sex) + (1 | qall) + (1 | del1p) + (1 | del17p) + (1 | trans_14_16) + (1 | trans_4_14)

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, fish_chip)
colnames(C) <- c("Stage", "Sex", "1q aberration", "1p del", "17p del", "t(14;16)", "t(4;14)")
rownames(C) <- colnames(C)  

# Plot correlation matrix
# between all pairs of variables

setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
png("correlation_matrix_noscore.png", res = 300, width = 2000, height = 2000)
plotCorrMatrix(C)
dev.off()

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,

varPart <- fitExtractVarPartModel(norm_data, form, fish_chip)

#Rturns an object that stores the variance fractions for each gene and each variable in the formula specified. These fractions can be accessed just like a data.frame:


# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
colnames(vp) <- c("17p del",  "1p del", "1q aber", "Sex", "Stage", "t(14;16)", "t(4;14)", "Residuals")
vp <- vp[,c("Stage", "t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")]

# Figure 1b
# violin plot of contribution of each variable to total variance
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
png("variance_explained_noscore.png", res = 300, width = 4000, height = 2000)
plotVarPart(vp)
dev.off()

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix/Model_testing")
save(varPart, file = "varPart_data_noscore.RData")