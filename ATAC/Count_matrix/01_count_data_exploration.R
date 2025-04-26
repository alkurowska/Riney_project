#####################
####  Prep DATA  ####
#####################

# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata") 
counts <- as.data.frame(Counts$counts) # 157817 x  248

# Load consensus peaks
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == consensus_peaks$GeneID)
# TRUE
rownames(counts) <- consensus_peaks$name # change rownames of the count data

# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
counts <- counts[,-which(colnames(counts) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]

# change colnames of the count data
# Remove everything after "_S"
colnames(counts) <- gsub("_S.*", "", colnames(counts))

# Remove SEX chromosomes
counts <- counts[!grepl("chrX|chrY", rownames(counts)),] # 152912  x  247
dim(counts)

# load metadata
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
fish <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(fish) <- fish[1,]
fish <- fish[-1,]

# Remove FISH data with NA values
fish <- fish[!is.na(fish$`1q amp`),] # 202
dim(fish)

counts <- counts[,colnames(counts)%in%fish$`Sample ATACseq`]
dim(counts) # 152912  x  190

fish <- fish[fish$`Sample ATACseq` %in% colnames(counts),]
dim(fish) # 190 x 17

rownames(fish) <- fish$`Sample ATACseq`


#####################
#### FRIP SCORES ####
#####################

# Get the data from mapping
stats <- Counts$stat 
rownames(stats) <- stats$Status
stats <- stats[,-1] # remove first column with status names
# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
stats <- stats[,-which(colnames(stats) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]
colnames(stats) <- gsub("_S.*", "", colnames(stats))
stats <- as.data.frame(t(stats))
stats$Reads <- NA
stats$ConsFRIP <- NA

# Calculate the FRIP score of the consensus peaks
for (i in 1:nrow(stats)) {
  stats$Reads[i] <- sum(stats[i,1:14])
  stats$ConsFRIP[i] <- stats$Assigned[i]/stats$Reads[i]
}

stats$ConsFRIP <- round(stats$ConsFRIP,digits = 2)
stats$sample <- rownames(stats)

# Keep only the samples to use
stats <- stats[stats$sample %in% rownames(fish),]

fish$FRIP <- stats[rownames(fish),]$ConsFRIP

########################
#### GSVA_low score ####
########################

# Add GSVA scores
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
fish <- fish[fish$`Sample RNAseq`%in%rownames(scores),]
dim(fish) #185 x 18
fish$gsva <- scores[fish$`Sample RNAseq`, "GSVA_low"]

counts <- counts[,colnames(counts)%in%rownames(fish)]
dim(counts) # 152912 x 185

table(fish$Stage)

# HC MGUS   MM  SMM 
#  6   14  136   29 

#######################
####     logCPM    ####
#######################

### create DGEList object
library(edgeR)
y <- DGEList(counts=counts)

## Filtering with the default parameters
stage <- factor(fish$Stage, levels = c("HC","MGUS", "SMM", "MM"))
keep <- filterByExpr(y, group=stage, min.count=10, min.total.count=15) 
table(keep)
# FALSE   TRUE 
# 10776 142136 

# keep peaks that pass the filter and recalculate the library size.
y <- y[keep, , keep.lib.sizes=FALSE]

# Transform data to counts per milion
## logCPM
cpm_data <- cpm(y)
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
write.table(cpm_data, "cpm_data.txt", sep="\t", dec=".", quote=FALSE)
log_cpm <- cpm(y, log = TRUE, prior.count = 1)
write.table(log_cpm, "log_cpm.txt", sep="\t", dec=".", quote=FALSE)


#### For models testing ####
# design model
fish_atac <- fish[colnames(y$counts),]
# Relevant variables for the design matrix
stage = as.factor(fish_atac$Stage)
design <- model.matrix(~ stage, data = fish_atac)
colnames(design) <- gsub("stage", "", colnames(design))
v <- voom(y, design, plot=F)
norm_data <- v$E

# Save the data
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing")
save(norm_data, v, consensus_peaks, fish_atac, file = "voom_to_models_testing.RData")
############################

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/")
save(y, fish_atac, consensus_peaks, file = "norm_to_voom.RData")

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
pca_atac <- prcomp(t(log_cpm))

# Plot PCs variance

library(factoextra)
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/")
png("PCA_variance.png", width = 800, height = 600)
fviz_screeplot(pca_atac, addlabels = TRUE, ylim = c(0, 30))
dev.off()


var_prcomp <- pca_atac$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

stage <- factor(as.character(fish_atac$Stage),levels=c("HC",
                                                                   "MGUS",
                                                                   "SMM",
                                                                   "MM"))
Stage <- stage
mds <- data.frame(PC1=pca_atac$x[,1],PC2=pca_atac$x[,2], Stage=Stage)


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

ggsave(plot = p3, dpi = 600, filename = "human_pca.jpg", width = 250, height = 200, units = "mm")


############################
####     BATCH EFFECT   ####
#### removeBatchEffect  ####
############################


# Load cytogenetic data for covariates of interest to be included in the model
# Substitute "neutral" for "" in the covariates
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

dim(fish_atac) 
fish_atac <- fish_atac[colnames(y$counts),]

# Relevant variables for the design matrix
stage = as.factor(fish_atac$Stage)
amp = as.factor(fish_atac$qall)
del1p = as.factor(fish_atac$del1p)
del17p = as.factor(fish_atac$del17p)
trans_14_16 = as.factor(fish_atac$trans_14_16)
trans_4_14 = as.factor(fish_atac$trans_4_14)
sex = as.factor(fish_atac$Sex)

FRIP <- as.numeric(fish_atac$FRIP)
gsva <- as.numeric(fish_atac$gsva) 



# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16, data = fish_atac)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"

# Remove batch effect
log_cpm_noFRIP <- removeBatchEffect(log_cpm, covariates = as.matrix(FRIP), design = design)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
write.table(log_cpm_noFRIP, "log_cpm_noFRIP.txt", sep="\t", dec=".", quote=FALSE)


 
# GSVA with  MIN value for MM
fish_atac[fish_atac$Stage == "MM", "gsva"] <- min(fish_atac[fish_atac$Stage == "MM", "gsva"])
gsva <- as.numeric(fish_atac$gsva)

# design model
design <- model.matrix(~ stage + sex + amp + del1p + del17p + trans_4_14 + trans_14_16, data = fish_atac)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("sex", "", colnames(design))
colnames(design) <- gsub("amp1q ", "q", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[9] <- "trans_4_14"
colnames(design)[10] <- "trans_14_16"
colnames(design)[1] <- "Intercept"
# Remove batch effect
log_cpm_noBatch <- removeBatchEffect(log_cpm, covariates = as.matrix(cbind(FRIP,gsva)), design = design)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
write.table(log_cpm_noBatch, "log_cpm_noBatch.txt", sep="\t", dec=".", quote=FALSE)


### STEP 3: UMAP + t-SNE
library(Rtsne)
library(umap)

# TSNE
atac_tsne <- Rtsne(t(log_cpm_noBatch), perplexity=20)

# UMAP
atac_umap <- umap(t(log_cpm_noBatch))

mplot <- data.frame(tSNE1=atac_tsne$Y[,1], tSNE2=atac_tsne$Y[,2], UMAP1=atac_umap$layout[,1], UMAP2=atac_umap$layout[,2], Stage=Stage)

p3 <- ggplot(mplot, aes(x = tSNE1, y = tSNE2, color = Stage)) +
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

p3 <- ggplot(mplot, aes(x = UMAP1, y = UMAP2, color = Stage)) +
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



#######################
####       PCA     ####
#######################

# PCA
pca_data <- prcomp(t(log_cpm_noBatch))

png("PCA_variance_regress.png", width = 800, height = 600)
fviz_screeplot(pca_data , addlabels = TRUE, ylim = c(0, 30))
dev.off()

# logCPM PCA
var_prcomp <- pca_data $sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)
library(RColorBrewer)

stage = as.factor(fish_atac$Stage)
mds <- data.frame(PC1=pca_data $x[,1],PC2=pca_data $x[,2], Stage=stage)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = -PC1, y = PC2, color = Stage)) +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_regress.jpg", width = 250, height = 200, units = "mm")



# PCA
pca_data <- prcomp(t(log_cpm_noFRIP))

library(factoextra)
png("PCA_variance_frip.png", width = 800, height = 600)
fviz_screeplot(pca_data , addlabels = TRUE, ylim = c(0, 30))
dev.off()

# logCPM PCA
var_prcomp <- pca_data $sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)
library(RColorBrewer)

stage = as.factor(fish_atac$Stage)
mds <- data.frame(PC1=pca_data $x[,1],PC2=pca_data $x[,2], Stage=stage)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = -PC1, y = -PC2, color = Stage)) +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_frip.jpg", width = 250, height = 200, units = "mm")
