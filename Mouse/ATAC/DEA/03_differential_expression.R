###########################################
####   RNA-Seq pipeline - SINGLE-END   ####
####   03_differential_expression.R    ####
####               LIMMA               ####
###########################################


# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("norm_to_voom.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
d1_norm <- d1_norm

#--------------------------------------------------------

### STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS >> VOOM+LIMMA

#Differential expression analysis
#require(edgeR)
combined <- as.factor(metadata$Combined)
batch <- as.factor(metadata$Batch)
# According to sex QC all unknown samples are female -> reassing 
metadata[metadata$Sex == "Unknown",]$Sex <- "female"
sex <- as.factor(metadata$Sex)

rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[colnames(d1_norm$counts),]
frip <- metadata$FRiP

design <- model.matrix(~ combined + sex + frip)
colnames(design) <- gsub("combined", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design)[1] <- "Intercept"

#voom transformation
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA")
png("voom.png", width = 800, height = 600)
v <- voom(d1_norm, design, plot=TRUE)
dev.off()

norm_data <- v$E
dim(norm_data) # 51685   49

# Save the data
save(norm_data, v, metadata, file = "voom_to_dea.RData")

# counts for samples with fish
atac_data <- norm_data[,rownames(metadata)] 
metadata <- metadata[colnames(atac_data),] 

# Differential expression analysis
fit <- lmFit(atac_data, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))


# # Check the difference between samples with 1q amp, 1q gain and 1q gainORamp for SMM and MM samples 
myContrMatrix = makeContrasts(
   CyclinD1_MM = CyclinD1_BIc_MM , #trasitiion
   CyclinD1_MGUS = CyclinD1_BIc_MGUS  , #trasitiion
   CyclinD1 = CyclinD1_BIc_MM - CyclinD1_BIc_MGUS , #trasitiion

    MIc_MM = MIc_MM  , #trasitiion
    MIc_MGUS = MIc_MGUS  , #trasitiion
    MIc = MIc_MM - MIc_MGUS , #trasitiion

    Mmset_MM = Mmset_BIc_MM  , #trasitiion
    Mmset_MGUS = Mmset_BIc_MGUS  , #trasitiion
    Mmset = Mmset_BIc_MM - Mmset_BIc_MGUS , #trasitiion

    Trp53_MM = Trp53_BIc_MM  , #trasitiion
    Trp53_MGUS = Trp53_BIc_MGUS  , #trasitiion
    Trp53 = Trp53_BIc_MM - Trp53_BIc_MGUS , #trasitiion
   
   levels= colnames(design)
  )

## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(atac_data), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #13991
sum(contrast_stats$adj.P.Val < 0.01) #5728

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

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0])  #4152
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #4152
#with a adj.pval<0.01 & |FC|>1.5 (logFC>0.58)

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA")
write.table(dea_results,"atac_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential genes per transition
length(which(dea_results$sig_CyclinD1_MM != 0)) # 2017
length(which(dea_results$sig_CyclinD1_MGUS != 0)) # 700
length(which(dea_results$sig_CyclinD1 != 0)) # 7

length(which(dea_results$sig_MIc_MM != 0)) # 1568
length(which(dea_results$sig_MIc_MGUS != 0)) # 1454
length(which(dea_results$sig_MIc != 0)) # 1

length(which(dea_results$sig_Mmset_MM != 0)) # 2122
length(which(dea_results$sig_Mmset_MGUS != 0)) # 421
length(which(dea_results$sig_Mmset != 0)) # 512

length(which(dea_results$sig_Trp53_MM != 0)) # 1868
length(which(dea_results$sig_Trp53_MGUS != 0)) # 977
length(which(dea_results$sig_Trp53 != 0)) # 1

# How many differential genes per transition are overlapping across models
# CyclinD1