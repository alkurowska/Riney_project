#########################################
# RNA-Seq pipeline - SINGLE-END         #
#                                       #
# script: 02_count_table_exploration.R  #
#########################################

# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")



####################### DATA PREPARATION #######################

### STEP 1: READ DATA 

# Library required
library(edgeR)
library(ggplot2)

## Count table
d <- read.table("count_table_NEW.txt", header = FALSE)
rownames(d) <- paste(d[,1],d[,2],d[,3],sep="_")
consensus_peaks <- d[,c(1:6)]
colnames(consensus_peaks) <- c("chr","start","end","name","score","strand")
counts <- d[,-c(1:6)]

# save consensus peaks
write.table(consensus_peaks, "consensus_peaks.txt", sep="\t", dec=".", quote=FALSE, row.names=T, col.names=TRUE)

#Name samples
lst <- c('249_pB2IC_Cp_Mieloma_S7.sort.rmdup.rmblackls.rmchr.bam',
         '2795_MmsetBIc_MGUS_S3.sort.rmdup.rmblackls.rmchr.bam',
         '2981_cyclinD1BIc_MGUS_S1.sort.rmdup.rmblackls.rmchr.bam',
         '2986_MmsetBIc_MGUS_S6.sort.rmdup.rmblackls.rmchr.bam',
         '3011_MmsetBIc_MGUS_S2.sort.rmdup.rmblackls.rmchr.bam',
         '3015_MmsetBIc_MGUS_S5.sort.rmdup.rmblackls.rmchr.bam',
         '3016_MmsetBIc_MGUS_S4.sort.rmdup.rmblackls.rmchr.bam',
         '3069_CyclinD1BIc_MM_S21.sort.rmdup.rmblackls.rmchr.bam',
         '3086_MmsetBIc_MM_S17.sort.rmdup.rmblackls.rmchr.bam',
         '3162_CyclinD1BIc_MM_S20.sort.rmdup.rmblackls.rmchr.bam',
         '3216_CyclinD1BIc_MGUS_S8.sort.rmdup.rmblackls.rmchr.bam',
         '3408_MmsetBIc_MM_S19.sort.rmdup.rmblackls.rmchr.bam',
         '3410_MmsetBIc_MM_S13.sort.rmdup.rmblackls.rmchr.bam',
         '3509_CyclinD1BIc_MGUS_S9.sort.rmdup.rmblackls.rmchr.bam',
        '3661_CyclinD1BIc_MM_S25.sort.rmdup.rmblackls.rmchr.bam',
         '3664_CyclinD1BIc_MGUS_S10.sort.rmdup.rmblackls.rmchr.bam',
         '3670_MmsetBIc_MM_S18.sort.rmdup.rmblackls.rmchr.bam',
         '3818_CyclinD1BIc_MGUS_S11.sort.rmdup.rmblackls.rmchr.bam',
         '3821_CyclinD1BIc_MGUS_S7.sort.rmdup.rmblackls.rmchr.bam',
         '4220_CyclinD1BIc_MGUS_S12.sort.rmdup.rmblackls.rmchr.bam',
         '5263_B2IC_Cp_Mieloma_S8.sort.rmdup.rmblackls.rmchr.bam',
         '5278_CyclinD1BIc_MM_S22.sort.rmdup.rmblackls.rmchr.bam',
         '6174_pB2IC_Cp_Mieloma_S9.sort.rmdup.rmblackls.rmchr.bam',
         '6434_B2IKC_Cp_Mieloma_S13.sort.rmdup.rmblackls.rmchr.bam',
         '6553_MmsetBIc_MM_S16.sort.rmdup.rmblackls.rmchr.bam',
         '6575_CyclinD1BIc_MM_S23.sort.rmdup.rmblackls.rmchr.bam',
         '6719_MmsetBIc_MM_S15.sort.rmdup.rmblackls.rmchr.bam',
         '6723_MmsetBIc_MM_S14.sort.rmdup.rmblackls.rmchr.bam',
         '6844_CyclinD1BIc_MM_S24.sort.rmdup.rmblackls.rmchr.bam',
         '8322_pB2IC_Cp_MGUS_S29.sort.rmdup.rmblackls.rmchr.bam',
         '8323_pB2IC_Cp_Mieloma_S31.sort.rmdup.rmblackls.rmchr.bam',
         '8326_pB2IC_Cp_MGUS_S32.sort.rmdup.rmblackls.rmchr.bam',
         '8327_pB2IC_Cp_Mieloma_S19.sort.rmdup.rmblackls.rmchr.bam',
         '8328_pB2IC_Cp_MGUS_S23.sort.rmdup.rmblackls.rmchr.bam',
         '8329_pB2IC_Cp_Mieloma_S11.sort.rmdup.rmblackls.rmchr.bam',
         '8330_pB2IC_Cp_MGUS_S24.sort.rmdup.rmblackls.rmchr.bam',
         '8331_pB2IC_Cp_MGUS_S18.sort.rmdup.rmblackls.rmchr.bam',
         '8332_pB2IC_Cp_Mieloma_S25.sort.rmdup.rmblackls.rmchr.bam',
         '8476_MIC_Mieloma_S15.sort.rmdup.rmblackls.rmchr.bam',
         '8582_MIC_Mieloma_S5.sort.rmdup.rmblackls.rmchr.bam',
         '8799_Maf_B2IKC_Mieloma_S26.sort.rmdup.rmblackls.rmchr.bam',
         '8936_MIC_MGUS_S4.sort.rmdup.rmblackls.rmchr.bam',
         '8982_MIC_MGUS_S3.sort.rmdup.rmblackls.rmchr.bam',
         '8985_MIC_MGUS_S1.sort.rmdup.rmblackls.rmchr.bam',
         '8991_MIC_MIeloma_S30.sort.rmdup.rmblackls.rmchr.bam',
         '8994_MIC_MIeloma_S28.sort.rmdup.rmblackls.rmchr.bam',
         '8995_MIC_MGUS_S2.sort.rmdup.rmblackls.rmchr.bam',
         '9312_Maf_MIC_Mieloma_S16.sort.rmdup.rmblackls.rmchr.bam',
         '9399_Maf_MIC_Mieloma_S17.sort.rmdup.rmblackls.rmchr.bam',
         'Pool6_YFPCg1_Control_S6.sort.rmdup.rmblackls.rmchr.bam',
         'Pool7_YFPCg1_Control_S10.sort.rmdup.rmblackls.rmchr.bam',
         'Pool9_YFP_Control_S22.sort.rmdup.rmblackls.rmchr.bam',
         'SC20_B2IC_Cp_Mieloma_S12.sort.rmdup.rmblackls.rmchr.bam',
         'SC30_MIC_Mieloma_S14.sort.rmdup.rmblackls.rmchr.bam',
         'SC40_MIC_Mieloma_S20.sort.rmdup.rmblackls.rmchr.bam',
         'SC60_Maf_MIC_Mieloma_S21.sort.rmdup.rmblackls.rmchr.bam',
         'SC80_MIC_MGUS_S27.sort.rmdup.rmblackls.rmchr.bam')
  

colnames(counts) <- lst
colnames(counts) <- gsub(".sort.rmdup.rmblackls.rmchr.bam","", colnames(counts))
dim(counts) #59089 x 57

#metadata
metadata <- read.table("Col_Data.csv",sep=",", header=TRUE, check.names=FALSE) #104 x 6
dim(metadata) #57 x 17

# -------------------------------

#KEEP ONLY 4 MODELS 

#remove samples from other models 
samples_excluded <- c("9399_Maf_MIC_Mieloma_S17",
                      "8799_Maf_B2IKC_Mieloma_S26",
                      "SC60_Maf_MIC_Mieloma_S21",
                      "6434_B2IKC_Cp_Mieloma_S13",
                      "9312_Maf_MIC_Mieloma_S16",
                      "SC20_B2IC_Cp_Mieloma_S12",
                      "5263_B2IC_Cp_Mieloma_S8")

counts <- counts[,!colnames(counts)%in%samples_excluded] #59089 x 50
dim(counts) 
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[colnames(counts),]
counts <- counts[,rownames(metadata)]

table(metadata$Sex)
# Remove SEX chromosomes
counts <- counts[!grepl("chrX|chrY", rownames(counts)),] # 59089  x  50
dim(counts)

### STEP 2: LIMMA VOOM NORMALIZATION 
require(edgeR)

##There are two strange samples: 
## SC40_MIC_Mieloma_S20

# Remove outliers
outliers <- c("SC40_MIC_Mieloma_S20")
counts2 <- counts[,!colnames(counts)%in%outliers] #59089 x 49
metadata <- metadata[!metadata$Sample_ID%in%outliers,] #49 x 6

metadata <- metadata[colnames(counts2),]

table(metadata$Model)
#Control CyclinD1_BIc          MIc    Mmset_BIc    Trp53_BIc 
#      3           12           10           13           11 
table(metadata$Stage)
# Control    MGUS      MM 
#       3      22      24 

#create DGEList object
d0 <- DGEList(counts=counts2)

#filter: remove rows that consistently have zero or very low counts
combined <- as.factor(metadata$Combined)
metadata[metadata$Batch == "1",]$Batch <- "batch_1"
metadata[metadata$Batch == "2",]$Batch <- "batch_2"
batch <- as.factor(metadata$Batch)
stage <- as.factor(metadata$Stage)
model <- as.factor(metadata$Model)

keep <- filterByExpr(d0, group=stage, min.count = 10, min.total.count=15) #51685
table(keep)
d1 <- d0[keep,,keep.lib.sizes=FALSE]

#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)
log_cpm <- cpm(d1_norm, log = TRUE, prior.count = 1)
cpm <- cpm(d1_norm, log = FALSE, prior.count = 1)

getwd()
saveRDS(d1$counts,"ATAC_raw_counts.RDS")
write.table(log_cpm, "log_cpm.txt", sep="\t", dec=".", quote=FALSE)
write.table(cpm, "cpm.txt", sep="\t", dec=".", quote=FALSE)

#voom transformation
design <- model.matrix(~ stage, data = metadata)
colnames(design) <- gsub("stage", "", colnames(design))
v <- voom(d1_norm, design, plot=F)
norm_data <- v$E

# Save the data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix/Model_testing")
save(norm_data, v, metadata, file = "voom_to_models_testing.RData")
############################

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix/")
save(d1_norm, metadata, file = "norm_to_voom.RData")

####################### EXPLORATORY ANALYSIS #######################
### STEP 3: PCA
library(ggplot2)
library(RColorBrewer)

# Labels colors
color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "SMM", "MM")

color_models <- c("#368716", "#CAB2DC", "#1F78B4",  "#FB9A99", "#CAF2B0")
names(color_models) <- levels(as.factor(metadata$Model))

color_batch <- c("#000000", "#D55E00", "#0072B2")
names(color_batch) <- levels(as.factor(metadata$Batch))

####################### EXPLORATORY ANALYSIS #######################
##Plot the PCA 
pca_rna <- prcomp(t(log_cpm)) 

library(factoextra)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix/")
png("PCA_variance.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 20))
dev.off()

var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

stage <- factor(as.character(metadata$Stage),levels=c("Control",
                                                                   "MGUS",
                                                                   "MM"))
Stage <- stage
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Stage=Stage)


p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = Stage)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #geom_text(label=rownames(mds), size = 2, check_overlap = F, vjust = -0.8, hjust = -0.1) + 
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca.jpg", width = 250, height = 200, units = "mm")

# Color by model
model <- factor(as.character(metadata$Model),levels=c("Control","Trp53_BIc", "MIc", "CyclinD1_BIc", "Mmset_BIc"))
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], PC3=pca_rna$x[,3],PC4=pca_rna$x[,4], 
                  stage=metadata$Stage, batch=metadata$Batch, model=model)

p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = model)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_models) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #geom_text(label=rownames(mds), size = 2, check_overlap = F, vjust = -0.8, hjust = -0.1) + 
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_model.jpg", width = 250, height = 200, units = "mm")


p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = batch)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_batch) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #geom_text(label=rownames(mds), size = 2, check_overlap = F, vjust = -0.8, hjust = -0.1) + 
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_batch.jpg", width = 250, height = 200, units = "mm")


# Check FRIP score 
frip <- metadata$FRiP

mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], PC3=pca_rna$x[,3],PC4=pca_rna$x[,4], 
                  stage=metadata$Stage, batch=metadata$Batch, model=model, frip =frip)

p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = frip)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_batch) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #geom_text(label=rownames(mds), size = 2, check_overlap = F, vjust = -0.8, hjust = -0.1) + 
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_frip.jpg", width = 250, height = 200, units = "mm")