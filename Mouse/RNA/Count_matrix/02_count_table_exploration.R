#########################################
# RNA-Seq pipeline - SINGLE-END         #
#                                       #
# script: 02_count_table_exploration.R  #
#########################################

# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")



####################### DATA PREPARATION #######################

### STEP 1: READ DATA 

#count matrix
load("RunsAll.rda")

#genes as row names
counts <- df 
rownames(counts) <- counts[, 1]
counts <- counts[,-1]

#remove the first 5 rows with the alignment information
counts <- counts[-c(1:5),] #55,419 x 104

#metadata
metadata <- read.table("Col_Data.csv",sep=",", header=TRUE, check.names=FALSE) #104 x 6

# -------------------------------

#KEEP ONLY 4 MODELS 

#remove samples from other models 
colnames(counts)
samples_excluded <- c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24",
                      "1108_B2IC_MGUS_S3", "1121_B2IC_MGUS_S4", "1148_B2IC_MGUS_S2",
                      "2358_KBIc_MM_S34", "3042_cMafB2IC_MM_S16", "3048_cMafB2IC_MM_S16",
                      "3253_B2IKC_MM_S1", "3313_cMafB2IC_MGUS_S8", "3498_cMafB2IC_MGUS_S8", "3502_cMafB2IC_MGUS_S9",
                      "3903_B2IKC_MGUS_S5", "3910_B2IKC_MGUS_S4", "3911_B2IC_MGUS_S3", "3940_B2IKC_MGUS_S5",            
                      "4294_KBIc_MM_S23", "5263_BIc_MM_S18", "5937_B2IC_MM_S14",           
                      "6430_KBIc_MM_S16", "6434_KBIc_MM_S14", "6538_B2IC_MM_S13",              
                      "8074_cMafB2IC_MM_S17", "8192_cMafB2IC_MM_S17", "8799_MafKBIc_MM_S1",
                      "9312_MafMIc_MM_S8", "9399_MafMIc_MM_S30", "SC20_Bic_MM_S6","SC60_MafMIc_MM_S11")

counts <- counts[,!colnames(counts)%in%samples_excluded] #55414 x 74
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[colnames(counts),]
counts <- counts[,metadata$Sample_ID]


### STEP 2: GENE ANNOTATION USING biomaRt

library(biomaRt)

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes_ids <- sub('\\.[0-9]*$', '', c(rownames(counts)))

biotypes = getBM(attributes = c("ensembl_gene_id", 'mgi_symbol','chromosome_name','start_position','end_position',
                                'strand','gene_biotype','percentage_gene_gc_content'), 
                 filters = "ensembl_gene_id", values=genes_ids, mart=mart) #54828 x 8 

#there are two mgi symbols for three ensembl gene id
which(duplicated(biotypes$ensembl_gene_id)) #[1] 24087 48622 52599
biotypes[24087,] 
biotypes[48622,]
biotypes[52599,]

biotypes2 <- biotypes[-which(biotypes$ensembl_gene_id=="ENSMUSG00000082414" & biotypes$mgi_symbol!="Gm21093"),]
biotypes2 <- biotypes2[-which(biotypes2$ensembl_gene_id=="ENSMUSG00000115016" & biotypes2$mgi_symbol!="Gm33906"),]
biotypes2 <- biotypes2[-which(biotypes2$ensembl_gene_id=="ENSMUSG00000119828" & biotypes2$mgi_symbol!="Gm53179"),]
#54825 x 8


sum(!(genes_ids %in% biotypes$ensembl_gene_id))
#there are 589 genes from counts table without annotation
#Check what type of ensembl ID are...
ids <- c(genes_ids[!genes_ids %in% biotypes2$ensembl_gene_id])

#remove those genes without annotation
rownames(counts) <- gsub("\\..*", "", rownames(counts))
counts_filtered <- counts[!rownames(counts)%in%ids,] #54,825 x 74

write.table(biotypes2,"gene_annotation.txt", sep="\t", dec=".", quote=FALSE, row.names=FALSE, col.names=TRUE)


### STEP 2: REMOVE Ig genes with chain IGH, IGK and IGL 
biotypes2 <- read.table("gene_annotation.txt", sep="\t", check.names = T)
colnames(biotypes2) <- biotypes2[1,]
biotypes2 <- biotypes2[-1,]

Igs <- biotypes2[grep("^Igh|Igk|Igl", biotypes2$mgi_symbol),]

# Check what kind of genes were retrieved 
levels(as.factor(Igs$gene_biotype))
# [1] "IG_C_gene"       "IG_D_gene"       "IG_D_pseudogene" "IG_J_gene"       "IG_V_gene"       "IG_V_pseudogene" "protein_coding"
# Check the protein coding genes
Igs_pro <- Igs[Igs$gene_biotype == "protein_coding",]
#.     ensembl_gene_id mgi_symbol chromosome_name start_position end_position strand   gene_biotype percentage_gene_gc_content
#1425  ENSMUSG00000013367     Iglon5               7       43122328     43139499     -1 protein_coding                      52.08
#4700  ENSMUSG00000024831    Ighmbp2              19        3309076      3333017     -1 protein_coding                      52.10
#20905 ENSMUSG00000075370      Igll1              16       16678535     16681849     -1 protein_coding                      47.63

#Remove Ighmbp2 and Igll1 from the list
Igs <- Igs[!Igs$mgi_symbol == "Ighmbp2",] # 395 x 8
Igs <- Igs[!Igs$mgi_symbol == "Igll1",] # 394 x 8

# add Jchain gene to the Ig list
jchain <- biotypes2[grep("^Jchain", biotypes2$mgi_symbol),]
Igs <- rbind(Igs, jchain) # 395 x 8 

write.table(Igs,"Ig_genes.txt", sep="\t", dec=".", quote=FALSE, row.names=FALSE, col.names=TRUE)


# Remove Immunoglobulins from count matrix 
counts2 <- counts_filtered[!rownames(counts_filtered)%in%Igs$ensembl_gene_id,] #54430 x 74
dim(counts2)

### STEP 2: LIMMA VOOM NORMALIZATION 
require(edgeR)

##There are two strange samples: 
## 3472_CD1B2IC_MGUS_S9
## 3872_MMSetB2IC_MGUS_S11

# Remove outliers
outliers <- c("3472_CD1B2IC_MGUS_S9", "3872_MMSetB2IC_MGUS_S11")
counts2 <- counts2[,!colnames(counts2)%in%outliers] #54430 x 72
metadata <- metadata[!metadata$Sample_ID%in%outliers,] #72 x 6

metadata <- metadata[colnames(counts2),]
#create DGEList object
d0 <- DGEList(counts=counts2)

#filter: remove rows that consistently have zero or very low counts
combined <- as.factor(metadata$Combined)
metadata[metadata$Batch == "1",]$Batch <- "batch_1"
metadata[metadata$Batch == "2",]$Batch <- "batch_2"
metadata[metadata$Batch == "3",]$Batch <- "batch_3"
batch <- as.factor(metadata$Batch)
stage <- as.factor(metadata$Stage)
model <- as.factor(metadata$Model)


keep <- filterByExpr(d0,group=stage, min.count = 1) #19423
table(keep)
d1 <- d0[keep,,keep.lib.sizes=FALSE]

#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)
log_cpm <- cpm(d1_norm, log = TRUE, prior.count = 1)
cpm <- cpm(d1_norm, log = FALSE, prior.count = 1)

getwd()
saveRDS(d1$counts,"RNA_raw_counts_noIGs.RDS")
write.table(log_cpm, "log_cpm.txt", sep="\t", dec=".", quote=FALSE)
write.table(cpm, "cpm.txt", sep="\t", dec=".", quote=FALSE)

#voom transformation
design <- model.matrix(~ stage, data = metadata)
colnames(design) <- gsub("stage", "", colnames(design))
v <- voom(d1_norm, design, plot=F)
norm_data <- v$E

# Save the data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix/Model_testing")
save(norm_data, v, biotypes2, metadata, file = "voom_to_models_testing.RData")
############################

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix/")
save(d1_norm, metadata, biotypes2, file = "norm_to_voom.RData")

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
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix/")
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

### STEP 6: BATCH CORRECTION
#-------------------------------------
## REMOVE BATCH EFFECT WITH LIMMA
batch <- as.factor(metadata$Batch)
stage <- as.factor(metadata$Stage)
model <- as.factor(metadata$Model)
combined <- as.factor(metadata$Combined)
design =  model.matrix(~ combined)
colnames(design) <- gsub("combined", "", colnames(design))

log_cpm_noBatch <- removeBatchEffect(log_cpm, batch = batch, design = design)
write.table(log_cpm_noBatch, "log_cpm_noBatch.txt", sep="\t", dec=".", quote=FALSE)

##Plot the PCA 
pca_data <- prcomp(t(log_cpm_noBatch))

png("PCA_variance_regress.png", width = 800, height = 600)
fviz_screeplot(pca_data , addlabels = TRUE, ylim = c(0, 25))
dev.off()

var_prcomp <- pca_data$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

stage <- factor(as.character(metadata$Stage),levels=c("Control",
                                                                   "MGUS",
                                                                   "MM"))
Stage <- stage
mds <- data.frame(PC1=pca_data$x[,1],PC2=pca_data$x[,2], Stage=Stage)


p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = Stage)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_regress.jpg", width = 250, height = 200, units = "mm")



mds <- data.frame(PC1=pca_data$x[,1],PC2=pca_data$x[,2], batch=batch, model=model, stage=stage)

p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = batch)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_batch) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_regress_batch.jpg", width = 250, height = 200, units = "mm")



p3 <- ggplot(mds, aes(x = -PC1, y = PC2, color = model)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_models) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_regress_model.jpg", width = 250, height = 200, units = "mm")


