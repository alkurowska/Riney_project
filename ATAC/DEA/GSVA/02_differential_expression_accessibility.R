###########################################
####         ATAC-Seq pipeline         ####
####  03_differential_accessibility.R  ####
####              HUMAN                ####
###########################################

set.seed(1234)

###########################


# Setting working directory
set.seed(1234)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("norm_to_voom.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
norm_data <- y
dim(norm_data)

# Load Cytogenetics 
# prepare the data for a design matrix
fish_atac <- fish_atac

### STEP 2: Multiple model analysis

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



#voom transformation
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
png("voom.png", width = 800, height = 600)
v <- voom(y, design, plot=TRUE)
dev.off()

norm_data <- v$E
dim(norm_data) # 142136 x 185

# Save the data

save(norm_data, v, consensus_peaks, fish_atac, file = "voom_to_dea.RData")


###################
####    DEA    ####
###################

atac_data <- norm_data
# Differential expression analysis
fit <- lmFit(atac_data, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))


# # Check the difference between each stage
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

contrast_stats <- topTable(fit.contrasts2, number=nrow(atac_data), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #25255
sum(contrast_stats$adj.P.Val < 0.01) #13087

dtest <- decideTests(fit.contrasts2)

contr.matrix <- myContrMatrix
all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(atac_data))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(atac_data))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(atac_data))

for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i,number=nrow(atac_data), sort.by="none")
  
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i] < 0.01),which(all_fch[,i] > 0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i] < 0.01),which(all_fch[,i] < -0.58))
  sig[sel2,i] <- -1 
}

colnames(all_padj) = paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) = paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) = paste("sig_",colnames(contr.matrix),sep="")

dea_results = data.frame(sig,all_padj,all_fch)
rownames(dea_results) = rownames(contrast_stats)
which(duplicated(rownames(dea_results)))

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0])  #22374
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #22374
#with a adj.pval<0.01 & |FC|>1.5 (logFC>0.58)

setwd("/ibex/user/kurowsaa/ATAC/DEA/GSVA")
write.table(dea_results,"atac_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential genes per transition
length(which(dea_results$sig_MM_SMM != 0)) # 0
length(which(dea_results$sig_SMM_MGUS != 0)) # 1
length(which(dea_results$sig_MGUS_HC != 0)) # 9540
length(which(dea_results$sig_MM_HC != 0)) # 20766
length(which(dea_results$sig_SMM_HC != 0)) # 16346
length(which(dea_results$sig_MM_MGUS != 0)) # 12



