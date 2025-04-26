###########################################
####     ATAC-Seq pipeline - DOUBLE    ####
####           Models_testing.R        ####
####              HUMAN                ####
###########################################


# Setting working directory
set.seed(1234)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("voom_to_models_testing.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
norm_data <- norm_data

# Load Cytogenetics 
fish_atac <- fish_atac
colnames(fish_atac)[4] <- "del17p"
colnames(fish_atac)[5] <- "del1p"
colnames(fish_atac)[6] <- "trans_4_14"
colnames(fish_atac)[7] <- "trans_14_16"
colnames(fish_atac)[11] <- "qall"

fish_atac$del17p <- gsub("neutral", "", fish_atac$del17p)
fish_atac$del1p <- gsub("neutral", "", fish_atac$del1p)
fish_atac$trans_4_14 <- gsub("neutral", "", fish_atac$trans_4_14)
fish_atac$trans_14_16 <- gsub("neutral", "", fish_atac$trans_14_16)
fish_atac$qall <- gsub("neutral", "", fish_atac$qall)

# Relevant covariates 
Stage = as.factor(fish_atac$Stage)
Sex = as.factor(fish_atac$Sex)
qall = as.factor(fish_atac$qall)
del1p = as.factor(fish_atac$del1p)
del17p = as.factor(fish_atac$del17p)
trans_14_16 = as.factor(fish_atac$trans_14_16)
trans_4_14 = as.factor(fish_atac$trans_4_14)


# Add GSVA scores
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
fish_atac <- fish_atac[fish_atac$`Sample RNAseq`%in%rownames(scores),]
dim(fish_atac) #185 x 18
fish_atac$gsva <- scores[fish_atac$`Sample RNAseq`, "GSVA_low"]

library("variancePartition")

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/fish_atac about each sample

# Stage together + 1q together
form <- ~ (1 | Stage) + (1 | Sex) + (1 | qall) + (1 | del1p) + (1 | del17p) + (1 | trans_14_16) + (1 | trans_4_14) + gsva

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, fish_atac)
colnames(C) <- c("Stage", "Sex", "1q aberration", "1p del", "17p del", "t(14;16)", "t(4;14)", "GSVA")
rownames(C) <- colnames(C)  

# Plot correlation matrix
# between all pairs of variables

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing")
png("correlation_matrix.png", res = 300, width = 2000, height = 2000)
plotCorrMatrix(C)
dev.off()

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,

varPart <- fitExtractVarPartModel(norm_data, form, fish_atac)

#Rturns an object that stores the variance fractions for each gene and each variable in the formula specified. These fractions can be accessed just like a data.frame:


# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)
colnames(vp) <- c("GSVA", "1q aber", "Stage", "17p del", "1p del", "Sex", "t(14;16)", "t(4;14)", "Residuals")
vp <- vp[,c("GSVA", "Stage","t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")]

# Figure 1b
# violin plot of contribution of each variable to total variance
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing")
png("variance_explained.png", res = 300, width = 4000, height = 2000)
plotVarPart(vp)
dev.off()

# Save data
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing")
save(varPart, file = "varPart_data.RData")