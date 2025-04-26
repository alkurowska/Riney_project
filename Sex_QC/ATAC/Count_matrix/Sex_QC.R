#####################
####  Prep DATA  ####
#####################

# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata") 
counts <- as.data.frame(Counts$counts) # 157817 x  248
colnames(counts)
# Load consensus peaks
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == consensus_peaks$GeneID)
# TRUE
rownames(counts) <- consensus_peaks$name # change rownames of the count data

# change colnames of the count data
# Remove everything after "_S"
colnames(counts) <- gsub("_S.*", "", colnames(counts))


# load metadata
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
fish <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(fish) <- fish[1,]
fish <- fish[-1,]

# Remove FISH data with NA values
fish <- fish[!is.na(fish$`1q amp`),] # 202
dim(fish)
colnames(fish)

atac_data <- counts

fish <- fish[fish$`Sample ATACseq` %in% colnames(atac_data),]
dim(fish) # 190 x 17

atac_data <- atac_data[,colnames(atac_data)%in%fish$`Sample ATACseq`]
dim(atac_data) # 157817  x  190

table(fish$Stage)
# HC MGUS   MM  SMM 
#  6   14  139   31 

rownames(fish) <- fish$`Sample ATACseq`
fish <- fish[colnames(atac_data),]

### STEP 2: LIMMA VOOM NORMALIZATION 

# Keep only peaks on Y chromosome
head(consensus_peaks)
toKeep <- consensus_peaks[consensus_peaks$Chr == "chrY",]$name # 723
cons_peaks <- consensus_peaks[consensus_peaks$Chr %in% c("chrY"),]

library(edgeR)
library(limma)
d0 <- DGEList(counts=atac_data, samples = fish, genes = consensus_peaks)

#filter: remove rows 
d1 <- d0[toKeep,,keep.lib.sizes=FALSE]

log_counts <- log2(d1$counts + 0.5)


####################### EXPLORATORY ANALYSIS #######################
### STEP 3: PCA
library(ggplot2)
library(RColorBrewer)


##Plot the PCA 
pca_rna <- prcomp(t(log_counts),center = T, scale = F)

# Plot PCs variance

library(factoextra)
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/ATAC/Count_matrix")
png("PCA_variance_sex_Y.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 80))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)
fish_rna <- fish


sex <- factor(as.character(fish$Sex),levels=c("female","male"))
Sex <- sex
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Sex = Sex)

# label only male samples 
mds$label <- ''
mds$label[mds$Sex == "female"] <- rownames(mds[mds$Sex == "female",])
library(ggrepel)
p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = mds$label), size = 2, show.legend = FALSE, max.overlaps = 40) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_Y_females.jpg", width = 250, height = 200, units = "mm")


# label only female samples 
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Sex = Sex)
mds$label <- ''
mds$label[mds$Sex == "male"] <- rownames(mds[mds$Sex == "male",])
library(ggrepel)
p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = mds$label), size = 2, show.legend = FALSE, max.overlaps = 40) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_Y_males.jpg", width = 250, height = 200, units = "mm")







# Keep only genes on X and Y chromosome
toKeep <- consensus_peaks[consensus_peaks$Chr %in% c("chrX", "chrY"),]$name # 4905

#filter: remove rows 
d1 <- d0[toKeep,,keep.lib.sizes=FALSE]

log_counts <- log2(d1$counts + 0.5)


####################### EXPLORATORY ANALYSIS #######################
### STEP 3: PCA
library(ggplot2)
library(RColorBrewer)


##Plot the PCA 
pca_rna <- prcomp(t(log_counts),center = T, scale = F)

# Plot PCs variance

library(factoextra)
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/ATAC/Count_matrix")
png("PCA_variance_sex_X_Y.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 80))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

sex <- factor(as.character(fish$Sex),levels=c("female","male"))
Sex <- sex
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], PC3=pca_rna$x[,3],Sex = Sex)

p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = rownames(mds)), size = 2, show.legend = FALSE) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_X_Y.jpg", width = 250, height = 200, units = "mm")




# label only male samples 
mds$label <- ''
mds$label[mds$Sex == "male"] <- rownames(mds[mds$Sex == "male",])
library(ggrepel)
p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = mds$label), size = 2, show.legend = FALSE, max.overlaps = 40) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_X_Y_males.jpg", width = 250, height = 200, units = "mm")


# label only female samples 
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Sex = Sex)
mds$label <- ''
mds$label[mds$Sex == "female"] <- rownames(mds[mds$Sex == "female",])
library(ggrepel)
p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = mds$label), size = 2, show.legend = FALSE, max.overlaps = 40) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_X_Y_females.jpg", width = 250, height = 200, units = "mm")
