###########################################
####     RNA-Seq pipeline - SINGLE     ####
####     03_differential_expression.R  ####
####              HUMAN                ####
###########################################


# Setting working directory
set.seed(1234)

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")

#Loading packages 
library(edgeR)
library(limma)

### STEP 1: LOAD DATA
load("norm_to_voom.RData")

#normalized data (Limma is using a linear model, data needs to be log transformed)
norm_data <- d1_norm
dim(norm_data)

# Load Cytogenetics 
# prepare the data for a design matrix
fish_rna <- fish_rna

### STEP 2: Multiple model analysis

colnames(fish_rna)[4] <- "del17p"
colnames(fish_rna)[5] <- "del1p"
colnames(fish_rna)[6] <- "trans_4_14"
colnames(fish_rna)[7] <- "trans_14_16"
colnames(fish_rna)[11] <- "qall"

fish_rna$qall <- gsub("neutral", "", fish_rna$qall)
fish_rna$del17p <- gsub("neutral", "", fish_rna$del17p)
fish_rna$del1p <- gsub("neutral", "", fish_rna$del1p)
fish_rna$trans_4_14 <- gsub("neutral", "", fish_rna$trans_4_14)
fish_rna$trans_14_16 <- gsub("neutral", "", fish_rna$trans_14_16)

# Relevant variables for the design matrix
stage = as.factor(fish_rna$Stage)
amp = as.factor(fish_rna$qall)
del1p = as.factor(fish_rna$del1p)
del17p = as.factor(fish_rna$del17p)
trans_14_16 = as.factor(fish_rna$trans_14_16)
trans_4_14 = as.factor(fish_rna$trans_4_14)
sex = as.factor(fish_rna$Sex)

# Add entropy and GSVA scores
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
table(rownames(fish_rna) == rownames(scores))
fish_rna$gsva <- scores[rownames(fish_rna), "GSVA_low"]
gsva <- as.numeric(fish_rna$gsva)

# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16 + gsva)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"



#voom transformation
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
png("voom.png", width = 800, height = 600)
v <- voom(d1_norm, design, plot=TRUE)
dev.off()

norm_data <- v$E
dim(norm_data) # 24256   197

# Save the data

save(norm_data, v, gene_anno, fish_rna, file = "voom_to_dea.RData")

# counts for samples with fish
rna_data <- norm_data[,rownames(fish_rna)] #197
fish_rna <- fish_rna[colnames(rna_data),] # 197

# Differential expression analysis
fit <- lmFit(rna_data, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))


# # Check the difference between samples with 1q amp, 1q gain and 1q gainORamp for SMM and MM samples 
myContrMatrix = makeContrasts(
   MM_SMM = MM - SMM, #trasitiion
   SMM_MGUS = SMM - MGUS, #trasitiion
   MM_MGUS = MM - MGUS, #trasitiion
   MGUS_HC = MGUS , #trasitiion
   MM_HC = MM , #trasitiion
   SMM_HC = SMM , #trasitiion
   levels= colnames(design)
  )

## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(rna_data), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #5514
sum(contrast_stats$adj.P.Val < 0.01) #2637

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

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0])  #5692
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #5692
#with a adj.pval<0.05 & |FC|>1.5 (logFC>0.58)

setwd("/ibex/user/kurowsaa/RNA/DEA/GSVA")
write.table(dea_results,"rna_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential genes per transition
length(which(dea_results$sig_MM_SMM != 0)) # 337
length(which(dea_results$sig_SMM_MGUS != 0)) # 36
length(which(dea_results$sig_MGUS_HC != 0)) # 1342
length(which(dea_results$sig_MM_HC != 0)) # 3247
length(which(dea_results$sig_SMM_HC != 0)) # 2664
length(which(dea_results$sig_MM_MGUS != 0)) # 2448



