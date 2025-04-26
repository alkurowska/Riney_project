###########################################
####     RNA-Seq pipeline - SINGLE     ####
####                                   ####
####              HUMAN                ####
###########################################

# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/RINEY/human/RNA_data")
library(edgeR)
library(limma)

####################### DATA PREPARATION #######################

### STEP 1: LOAD RAW DATA without Igs with healthy young donors
rna_data <- readRDS("RNA_Counts_NoFilter.RDS") #60252  x  216
dim(rna_data)

# Gene annotation 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/QC")
gene_anno <- readRDS("Gene_anno_noIGs.RDS") #60217 x 12
table(gene_anno$gene_id%in%rownames(rna_data))

# Filter rna data
rna_data <- rna_data[gene_anno$gene_id,]
dim(rna_data) # 60217 x 216

# load metadata
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
fish <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(fish) <- fish[1,]
fish <- fish[-1,]

# Remove FISH data with NA values
fish <- fish[!is.na(fish$`1q amp`),] # 202
dim(fish)

rna_data <- rna_data[,colnames(rna_data)%in%fish$`Sample RNAseq`]
dim(rna_data) # 60252 x 197

fish <- fish[fish$`Sample RNAseq` %in% colnames(rna_data),]
dim(fish) # 197 x 17

rownames(fish) <- fish$`Sample RNAseq`

table(fish$Stage)
# HC MGUS   MM  SMM 
#  7   15  145   30 

### STEP 2: LIMMA VOOM NORMALIZATION 

#create DGEList object
gene_anno <- gene_anno[gene_anno$gene_id %in% rownames(rna_data),] 
dim(gene_anno) # 60217 x 12

d0 <- DGEList(counts=rna_data, samples = fish, genes = gene_anno)

#filter: remove rows that consistently have zero or very low counts
stage <- as.factor(fish$Stage)

keep <- filterByExpr(d0, group = stage)
table(keep)
#keep
#FALSE  TRUE 
#35961 24256  
d1 <- d0[keep,,keep.lib.sizes=FALSE]


#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)
log_cpm <- cpm(d1_norm, log = TRUE, prior.count = 1)
cpm <- cpm(d1_norm, log = FALSE, prior.count = 1)

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
saveRDS(d1$counts,"RNA_raw_counts_noIGs.RDS")
write.table(log_cpm, "log_cpm.txt", sep="\t", dec=".", quote=FALSE)
write.table(cpm, "cpm.txt", sep="\t", dec=".", quote=FALSE)


#### For models testing ####
# design model
fish_rna <- fish[colnames(d1_norm$counts),]
# Relevant variables for the design matrix
stage = as.factor(fish_rna$Stage)
design <- model.matrix(~ stage, data = fish_rna)
colnames(design) <- gsub("stage", "", colnames(design))
v <- voom(d1_norm, design, plot=F)
norm_data <- v$E

# Save the data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
save(norm_data, v, gene_anno, fish_rna, file = "voom_to_models_testing.RData")
############################

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
save(d1_norm, fish_rna, gene_anno, file = "norm_to_voom.RData")

####################### EXPLORATORY ANALYSIS #######################
### STEP 3: PCA
library(ggplot2)
library(RColorBrewer)

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))
colors <- colors[-1]

color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("HC", "MGUS", "SMM", "MM")

##Plot the PCA 
pca_rna <- prcomp(t(log_cpm))

# Plot PCs variance

library(factoextra)
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
png("PCA_variance.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 10))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

stage <- factor(as.character(fish_rna$Stage),levels=c("HC",
                                                                   "MGUS",
                                                                   "SMM",
                                                                   "MM"))
Stage <- stage
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Stage=Stage)


p3 <- ggplot(mds, aes(x = PC1, y = -PC2, color = Stage)) +
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

ggsave(plot = p3, dpi = 600, filename = "human_pca.jpg", width = 250, height = 200, units = "mm")

############################
####     BATCH EFFECT   ####
#### removeBatchEffect  ####
############################

# Cytogenetic data for covariates of interest to be included in the model
# Substitute "neutral" for "" in the covariates

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

# Add GSVA scores
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
table(rownames(fish_rna) == rownames(scores))
fish_rna$entropy <- scores[rownames(fish_rna), "Igh_entropy"]
fish_rna$gsva <- scores[rownames(fish_rna), "GSVA_low"]
gsva <- as.numeric(fish_rna$gsva)

 
# GSVA with MIN value for MM
fish_rna[fish_rna$Stage == "MM", "gsva"] <- min(fish_rna[fish_rna$Stage == "MM", "gsva"])
gsva <- as.numeric(fish_rna$gsva)

# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"

saveRDS(design,"design_matrix.RDS")
# Remove batch effect
log_cpm_noBatch <- removeBatchEffect(log_cpm, covariates = gsva, design = design)
write.table(log_cpm_noBatch, "log_cpm_noBatch.txt", sep="\t", dec=".", quote=FALSE)

#######################
####       PCA     ####
#######################

# PCA
pca_data <- prcomp(t(log_cpm_noBatch))

png("PCA_variance_regress.png", width = 800, height = 600)
fviz_screeplot(pca_data , addlabels = TRUE, ylim = c(0, 15))
dev.off()

var_prcomp <- pca_data$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

stage <- factor(as.character(fish_rna$Stage),levels=c("HC",
                                                                   "MGUS",
                                                                   "SMM",
                                                                   "MM"))
Stage <- stage
mds <- data.frame(PC1=pca_data$x[,1],PC2=pca_data$x[,2], Stage=Stage)


p3 <- ggplot(mds, aes(x = PC1, y = -PC2, color = Stage)) +
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

ggsave(plot = p3, dpi = 600, filename = "human_pca_regress.jpg", width = 250, height = 200, units = "mm")



### STEP 3: UMAP + t-SNE
library(Rtsne)
library(umap)

# TSNE
rna_tsne <- Rtsne(t(log_cpm_noBatch), perplexity=20)

# UMAP
rna_umap <- umap(t(log_cpm_noBatch))

mplot <- data.frame(tSNE1=rna_tsne$Y[,1], tSNE2=rna_tsne$Y[,2], UMAP1=rna_umap$layout[,1], UMAP2=rna_umap$layout[,2], Stage=Stage)

p3 <- ggplot(mplot, aes(x = -tSNE1, y = -tSNE2, color = Stage)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(title="t-SNE - log CPM")
  # labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
  #   x = paste("t-SNE1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("t-SNE2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_tsne.jpg", width = 250, height = 200, units = "mm")

p3 <- ggplot(mplot, aes(x = -UMAP1, y = UMAP2, color = Stage)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(title="UMAP - log CPM")
  # labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
  #   x = paste("t-SNE1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("t-SNE2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_umap.jpg", width = 250, height = 200, units = "mm")
