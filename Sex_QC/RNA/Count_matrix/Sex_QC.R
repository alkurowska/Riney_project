###########################################
####     RNA-Seq pipeline - SINGLE     ####
####                                   ####
####              HUMAN                ####
###########################################



# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/RINEY/human/Stratification/")
library(edgeR)
library(limma)

####################### DATA PREPARATION #######################

### STEP 1: LOAD RAW DATA without Igs with healthy young donors
rna_data <- readRDS("RNA_Counts_NoFilter.RDS") #60252  x  216
dim(rna_data)

# Gene annotation 
gene_anno <- readRDS("Complete_GTF.RDS") #60666 x 12

# load metadata
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
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

fish <- fish[colnames(rna_data),]

### STEP 2: LIMMA VOOM NORMALIZATION 

#create DGEList object
gene_anno <- gene_anno[gene_anno$gene_id %in% rownames(rna_data),] # remove Igs
dim(gene_anno) # 60252 x 12

# Keep only genes on Y chromosome
toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #521

d0 <- DGEList(counts=rna_data, samples = fish, genes = gene_anno)

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
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
png("PCA_variance_sex_Y.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 80))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)
fish_rna <- fish


sex <- factor(as.character(fish_rna$Sex),levels=c("female","male", "unknown"))
Sex <- sex
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], Sex = Sex)

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

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_Y_males.jpg", width = 250, height = 200, units = "mm")


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

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_Y_females.jpg", width = 250, height = 200, units = "mm")







# Keep only genes on X and Y chromosome
toKeep <- gene_anno[gene_anno$V1 %in% c("Y", "X"),]$gene_id # #2945

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
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
png("PCA_variance_sex_X_Y.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 80))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

sex <- factor(as.character(fish_rna$Sex),levels=c("female","male", "unknown"))
Sex <- sex
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], PC3=pca_rna$x[,3],Sex = Sex)

library(ggrepel)
p3 <- ggplot(mds, aes(x = PC2, y = PC3, color = Sex)) +
  #geom_point(size = 6.5, show.legend = TRUE, shape = 19, color = "black") +
  #scale_color_manual(values = color_sex) + 
  geom_point(size = 2, show.legend = TRUE, shape = 19) +
  # add point labels names from mds rownames
  geom_text_repel(aes(label = rownames(mds)), size = 2, show.legend = FALSE) +

  #directlabels::geom_dl(aes(label = Stage), method = list("smart.grid", cex=1.5))+
  #ggforce::geom_mark_ellipse(aes(group = Stage, label = Stage)) + 
  labs(#title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",
    x = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""), y = paste("PC3 (",round(pcvar[3,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20))

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_X_Y.jpg", width = 250, height = 200, units = "mm")

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

ggsave(plot = p3, dpi = 600, filename = "human_pca_sex_X_Y_pc1pc2.jpg", width = 250, height = 200, units = "mm")




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




# BOX plot gene expression of chrY genes for the female samples, htat cluster close to male
# Keep only genes on Y chromosome
toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #521

female <- c("527_5", "300_03", "348_23", "344_7", "526_1", "543_3","790_1")
male <- fish$`Sample RNAseq`[fish$Sex == "male"]

samples <- c(female,male)

d1 <- d0[toKeep,,keep.lib.sizes=FALSE]
log_counts <- d1$counts

log_counts <- log_counts[,colnames(log_counts)%in%samples]

library(ggplot2)
library(RColorBrewer)
# Ensure `log_counts` is in a data frame format if not already
log_counts_df <- as.data.frame(log_counts)
log_counts_df$Gene <- rownames(log_counts_df)  # Add gene names if they're rownames
log_counts_long <- tidyr::pivot_longer(log_counts_df, cols = -Gene, names_to = "Sample", values_to = "Expression")

fish_plot <- fish[samples,]
colnames(fish_plot)[2] <- "Sample"
colnames(fish_plot)[1] <- "name"
# Merge sex information
log_counts_long <- merge(log_counts_long, fish_plot, by = "Sample")
log_counts_long$Sex <- factor(log_counts_long$Sex, levels = c("female", "male"))

# Order the Sample factor by Sex, then alphabetically within each Sex
log_counts_long$Sample <- factor(
  log_counts_long$Sample,
  levels = unique(log_counts_long$Sample[order(log_counts_long$Sex, log_counts_long$Sample)])
)


# Plot the boxplot
p3 <- ggplot(data = log_counts_long, aes(x = Sample, y = Expression, fill = Sex)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression of chrY Genes", x = "Samples", y = "Log2(Counts + 0.5)") +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  scale_fill_manual(values = c("male" = "blue", "female" = "pink"))  # Customize colors if needed

  
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
ggsave(plot = p3, dpi = 600, filename = "human_boxplot_chrY_genes.jpg", width = 1000, height = 200, units = "mm")



# Repeat for all of the chrX genes
toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #2424

male <- c("538_8", "812_05", "527_4", "348_28", "CUN_872774", "820_5")
rest <- fish$`Sample RNAseq`[fish$Sex == "male"]

samples <- unique(c(male,rest))
toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #521

d1 <- d0[toKeep,,keep.lib.sizes=FALSE]
log_counts <- log2(d1$counts + 0.5)

log_counts <- log_counts[,colnames(log_counts)%in%samples]

library(ggplot2)
library(RColorBrewer)
# Ensure `log_counts` is in a data frame format if not already
log_counts_df <- as.data.frame(log_counts)
log_counts_df$Gene <- rownames(log_counts_df)  # Add gene names if they're rownames
log_counts_long <- tidyr::pivot_longer(log_counts_df, cols = -Gene, names_to = "Sample", values_to = "Expression")

fish_plot <- fish[samples,]
colnames(fish_plot)[2] <- "Sample"
colnames(fish_plot)[1] <- "name"

# for this  samples male <- c("538_8", "812_05", "527_4", "348_28", "CUN_872774", "820_05") change Sex into "maybe female"
fish_plot[male,]$Sex <- "maybe female"
# remove NA
fish_plot <- fish_plot[!is.na(fish_plot$Stage),]

# Merge sex information
log_counts_long <- merge(log_counts_long, fish_plot, by = "Sample")
log_counts_long$Sex <- factor(log_counts_long$Sex, levels = c("maybe female", "male"))

# Order the Sample factor by Sex, then alphabetically within each Sex
log_counts_long$Sample <- factor(
  log_counts_long$Sample,
  levels = unique(log_counts_long$Sample[order(log_counts_long$Sex, log_counts_long$Sample)])
)


# Plot the boxplot
p3 <- ggplot(data = log_counts_long, aes(x = Sample, y = Expression, fill = Sex)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression of chrY Genes", x = "Samples", y = "Log2(Counts + 0.5)") +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  scale_fill_manual(values = c("male" = "blue", "maybe female" = "pink"))  # Customize colors if needed

  
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
ggsave(plot = p3, dpi = 600, filename = "human_boxplot_chrY_genes_males.jpg", width = 1000, height = 200, units = "mm")







# Repeat for all of the chrX genes
toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #2424

samples <- fish$`Sample RNAseq`[fish$Sex == "female"]

toKeep <- gene_anno[gene_anno$V1 == "Y",]$gene_id # #521

d1 <- d0[toKeep,,keep.lib.sizes=FALSE]
log_counts <- log2(d1$counts + 0.5)

log_counts <- log_counts[,colnames(log_counts)%in%samples]

library(ggplot2)
library(RColorBrewer)
# Ensure `log_counts` is in a data frame format if not already
log_counts_df <- as.data.frame(log_counts)
log_counts_df$Gene <- rownames(log_counts_df)  # Add gene names if they're rownames
log_counts_long <- tidyr::pivot_longer(log_counts_df, cols = -Gene, names_to = "Sample", values_to = "Expression")

fish_plot <- fish[samples,]
colnames(fish_plot)[2] <- "Sample"
colnames(fish_plot)[1] <- "name"

# remove NA
fish_plot <- fish_plot[!is.na(fish_plot$Stage),]

# Merge sex information
log_counts_long <- merge(log_counts_long, fish_plot, by = "Sample")
log_counts_long$Sex <- factor(log_counts_long$Sex, levels = c("female", "male"))

# Order the Sample factor by Sex, then alphabetically within each Sex
log_counts_long$Sample <- factor(
  log_counts_long$Sample,
  levels = unique(log_counts_long$Sample[order(log_counts_long$Sex, log_counts_long$Sample)])
)


# Plot the boxplot
p3 <- ggplot(data = log_counts_long, aes(x = Sample, y = Expression, fill = Sex)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression of chrY Genes", x = "Samples", y = "Log2(Counts + 0.5)") +
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  scale_fill_manual(values = c("male" = "blue", "female" = "pink"))  # Customize colors if needed

  
setwd("/ibex/user/kurowsaa/RINEY/Riney_pipeline/BULK/RNA/Count_matrix")
ggsave(plot = p3, dpi = 600, filename = "human_boxplot_chrY_genes_females.jpg", width = 1000, height = 200, units = "mm")


