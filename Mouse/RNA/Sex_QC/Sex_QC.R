###########################################
####     RNA-Seq pipeline - SINGLE     ####
####                                   ####
####              HUMAN                ####
###########################################



# Setting working directory
set.seed(123)
library(edgeR)
library(limma)

####################### DATA PREPARATION #######################

### STEP 1: LOAD RAW DATA without Igs 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")
rna_data <- readRDS("RNA_raw_counts_noIGs.RDS") # 19423 x 72
dim(rna_data)

# Gene annotation 
gene_anno <- read.table("gene_annotation.txt", sep = "\t", header = T)
head(gene_anno)

# load metadata
#metadata
metadata <- read.table("Col_Data.csv",sep=",", header=TRUE, check.names=FALSE) #104 x 6
metadata <- metadata[metadata$"Sample_ID" %in% colnames(rna_data),] # 72 x 6
rownames(metadata) <- metadata$Sample_ID

table(metadata$Stage)
# HC MGUS   MM 
#  6   32   34
table(metadata$Model)
# Control CyclinD1_BIc          MIc    Mmset_BIc    Trp53_BIc 
#       6           18           14           17           17 
table(metadata[,c("Stage", "Model")])
#           Control CyclinD1_BIc MIc Mmset_BIc Trp53_BIc
#   Control       6            0   0         0         0
#   MGUS          0            8   6         9         9
#   MM            0           10   8         8         8



### STEP 2: LIMMA VOOM NORMALIZATION 

#create DGEList object
gene_anno <- gene_anno[gene_anno$ensembl_gene_id %in% rownames(rna_data),] # remove Igs
dim(gene_anno) # 19423 x 8

# Keep only genes on X and Y chromosome
toKeep <- gene_anno[gene_anno$chromosome_name %in% c("Y", "X"),]$ensembl_gene_id # 635

d0 <- DGEList(counts=rna_data)

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
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Sex_QC")
png("PCA_variance_sex_X_Y.png", width = 800, height = 600)
fviz_screeplot(pca_rna, addlabels = TRUE, ylim = c(0, 40))
dev.off()


var_prcomp <- pca_rna$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

library(ggplot2)

# Check sample order
table(rownames(metadata) == colnames(rna_data))
# reorder 
metadata <- metadata[colnames(rna_data),]

sex <- factor(as.character(metadata$Sex),levels=c("female","male", "mix", "Unknown"))
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

ggsave(plot = p3, dpi = 600, filename = "mouse_pca_sex_X_Y.jpg", width = 250, height = 200, units = "mm")


### >>>>> All Unknown are female ?
