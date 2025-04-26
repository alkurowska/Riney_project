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

rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[colnames(d1_norm$counts),]

design <- model.matrix(~ combined + sex + batch)
colnames(design) <- gsub("combined", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("batchbatch", "batch", colnames(design))
colnames(design)[1] <- "Intercept"

#voom transformation
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
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

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0])  #5978
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #5978
#with a adj.pval<0.01 & |FC|>1.5 (logFC>0.58)

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
write.table(dea_results,"rna_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential genes per transition
length(which(dea_results$sig_CyclinD1_MM != 0)) # 3813
length(which(dea_results$sig_CyclinD1_MGUS != 0)) # 623
length(which(dea_results$sig_CyclinD1 != 0)) # 168

length(which(dea_results$sig_MIc_MM != 0)) # 2093
length(which(dea_results$sig_MIc_MGUS != 0)) # 871
length(which(dea_results$sig_MIc != 0)) # 0

length(which(dea_results$sig_Mmset_MM != 0)) # 3102
length(which(dea_results$sig_Mmset_MGUS != 0)) # 88
length(which(dea_results$sig_Mmset != 0)) # 793

length(which(dea_results$sig_Trp53_MM != 0)) # 3705
length(which(dea_results$sig_Trp53_MGUS != 0)) # 415
length(which(dea_results$sig_Trp53 != 0)) # 397

# How many differential genes per transition are overlapping across models
# CyclinD1
length(which(dea_results$sig_CyclinD1_MM == 1 & dea_results$sig_CyclinD1_MGUS == 1)) # 278
length(which(dea_results$sig_CyclinD1_MM == -1 & dea_results$sig_CyclinD1_MGUS == -1)) # 311

length(which(dea_results$sig_MIc_MM == 1 & dea_results$sig_MIc_MGUS == 1)) # 348
length(which(dea_results$sig_MIc_MM == -1 & dea_results$sig_MIc_MGUS == -1)) # 386

length(which(dea_results$sig_Mmset_MM == 1 & dea_results$sig_Mmset_MGUS == 1)) # 54
length(which(dea_results$sig_Mmset_MM == -1 & dea_results$sig_Mmset_MGUS == -1)) # 33

length(which(dea_results$sig_Trp53_MM == 1 & dea_results$sig_Trp53_MGUS == 1)) # 144
length(which(dea_results$sig_Trp53_MM == -1 & dea_results$sig_Trp53_MGUS == -1)) # 241

