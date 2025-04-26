###########################################
####   RNA-Seq pipeline - SINGLE-END   ####
####   03_differential_expression.R    ####
####               LIMMA               ####
###########################################


# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")

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
stage <- as.factor(metadata$Stage)
model <- as.factor(metadata$Model)

rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[colnames(d1_norm$counts),]

design <- model.matrix(~ stage + model + sex + batch)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("model", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("batchbatch", "batch", colnames(design))
colnames(design)[1] <- "Intercept"

#voom transformation
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/ALL_MODELS")
png("voom.png", width = 800, height = 600)
v <- voom(d1_norm, design, plot=TRUE)
dev.off()

norm_data <- v$E
dim(norm_data) # 19423   72

# Save the data
save(norm_data, v, biotypes2, metadata, file = "voom_to_dea.RData")

# counts for samples with fish
rna_data <- norm_data[,rownames(metadata)] 
metadata <- metadata[colnames(rna_data),] 

# Differential expression analysis
fit <- lmFit(rna_data, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))


# # Check the difference between samples with 1q amp, 1q gain and 1q gainORamp for SMM and MM samples 
myContrMatrix = makeContrasts(
   MM_MGUS = MM - MGUS, #trasitiion
   MGUS_HC = MGUS , #trasitiion
   MM_HC = MM , #trasitiion
   levels= colnames(design)
  )

## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(rna_data), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #7889
sum(contrast_stats$adj.P.Val < 0.01) #4768

dtest <- decideTests(fit.contrasts2)

contr.matrix <- myContrMatrix
all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(rna_data))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(rna_data))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(rna_data))

for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i,number=nrow(rna_data), sort.by="none")
  
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] > 0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] < -0.58))
  sig[sel2,i] <- -1 
}

colnames(all_padj) = paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) = paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) = paste("sig_",colnames(contr.matrix),sep="")

dea_results = data.frame(sig,all_padj,all_fch)
rownames(dea_results) = rownames(contrast_stats)
which(duplicated(rownames(dea_results)))

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0])  #7419
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #7419
#with a adj.pval<0.05 & |FC|>1.5 (logFC>0.58)

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/ALL_MODELS")
write.table(dea_results,"rna_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential genes per transition
length(which(dea_results$sig_MGUS_HC != 0)) # 2210
length(which(dea_results$sig_MM_HC != 0)) # 6527
length(which(dea_results$sig_MM_MGUS != 0)) # 2848

# How many differential genes per transition are overlapping 
length(which(dea_results$sig_MM_HC == 1 & dea_results$sig_MGUS_HC == 1)) # 579
length(which(dea_results$sig_MM_HC == -1 & dea_results$sig_MGUS_HC == -1)) # 1542

length(which(dea_results$sig_MM_HC == 1 & dea_results$sig_MM_MGUS == 1)) # 708
length(which(dea_results$sig_MM_HC == -1 & dea_results$sig_MM_MGUS == -1)) # 1320
