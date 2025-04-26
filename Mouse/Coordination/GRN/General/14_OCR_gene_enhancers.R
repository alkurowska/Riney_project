#Log cpm - tmm  - ATAC data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")
atac_matrix <- read.table("log_cpm.txt", header=T, sep="\t")  
colnames(atac_matrix) <- gsub("X", "", colnames(atac_matrix))
dim(atac_matrix) # 51685    49

#Log cpm - tmm  - RNA data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")
rna_matrix <- read.table("log_cpm.txt", header=TRUE)
colnames(rna_matrix) <- gsub("X", "", colnames(rna_matrix))
dim(rna_matrix) # 19423   72

# MATCH ATAC & RNA names
# read metadata
setwd("/ibex/user/kurowsaa/Riney_project/Mouse")
metadata <- read.table("metadata.csv", sep=",", header=TRUE)
head(metadata)

# Assign RNA names to ATAC
for(i in 1:ncol(atac_matrix)){
  atac_name <- colnames(atac_matrix)[i]
  new_name <- metadata[which(metadata$ATAC_ID==atac_name),]$Sample_ID
  colnames(atac_matrix)[i] <- new_name
}

# Paired samples between ATAC and RNA
common <- intersect(colnames(atac_matrix), colnames(rna_matrix))
atac_common <- atac_matrix[,common]
rna_common <- rna_matrix[,common]

dim(atac_common) # 51685    47
dim(rna_common) # 19423   47

table(colnames(atac_common) == colnames(rna_common)) # TRUE


#OCR genomic annotation
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Peaks")
peak_annotation <- readRDS("peak_annotation.RDS")
ocr_anno <- peak_annotation
rownames(ocr_anno) <-paste0(ocr_anno$seqnames, "_", ocr_anno$start, "_", ocr_anno$end)

ocr_anno <- ocr_anno[rownames(atac_common),]

# Load gene annotation
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")
gene_anno <- read.table("gene_annotation.txt", sep="\t", header=TRUE)
gene_anno <- gene_anno[gene_anno$ensembl_gene_id%in%rownames(rna_common),]
rownames(gene_anno) <- gene_anno$ensembl_gene_id
gene_anno <- gene_anno[rownames(rna_common),]

table(rownames(gene_anno) == rownames(rna_common) )# TRUE


#-------------------------------------

rna_common <- as.matrix(rna_common)
atac_common <- as.matrix(atac_common)



##>>>>>>>>>>>>>>>>>>>>> Correlation

# ENHANCERS 
anno_enhancer <- ocr_anno[which(ocr_anno$annotation2=="Distal Intergenic" | ocr_anno$annotation2=="Intron" ),]
dim(anno_enhancer) # 34200   20
# remove NA ENSEMBL
#anno_enhancer <- anno_enhancer[!is.na(anno_enhancer$ENSEMBL),]
#dim(anno_enhancer) 
atac_enhancer <- atac_common[rownames(anno_enhancer),] #intronic/distal intergenic peaks info
dim(atac_enhancer) # 34200   47


atac_postion <- strsplit(rownames(anno_enhancer),"_")
atac_start <- as.numeric(unlist(lapply(atac_postion,FUN=function(x){return(x[2])})))
atac_end <- as.numeric(unlist(lapply(atac_postion,FUN=function(x){return(x[3])})))
atac_chr <- gsub("chr","",unlist(lapply(atac_postion,FUN=function(x){return(x[1])})))
names(atac_start) <- rownames(atac_enhancer)
names(atac_end) <- rownames(atac_enhancer)

#Get the TSS and define the window
summary(rownames(gene_anno) == rownames(rna_common))
colnames(gene_anno)
gene_anno$TSS <- rep(NA, nrow(gene_anno))
for(i in 1:nrow(gene_anno)){
  if(gene_anno$strand[i]=="-1"){ gene_anno$TSS[i] <- gene_anno$end_position[i]} 
  if(gene_anno$strand[i]=="1"){ gene_anno$TSS[i] <- gene_anno$start_position[i]}
}
gene_anno$TSS <- as.numeric(gene_anno$TSS)

##Define the window of 500kb
window_up <- as.numeric(gene_anno$TSS)-500000
window_dw <- as.numeric(gene_anno$TSS)+500000


#--------------------------------------------------------------------#
#*********************    For all ENHANCERS   ***********************#
#--------------------------------------------------------------------#


#--------------------------------------------------------------------------#
#*********************          ALL SAMPLES         ***********************#
#--------------------------------------------------------------------------#
# Keep all samples

atac <- atac_common
rna <- rna_common

peak_gene_association <- list()

for(i in 1:nrow(rna)){ #for each gene i
# select peaks with that gene i
  sel <- which(atac_chr==gene_anno$chromosome_name[i] & atac_start >= window_up[i] & atac_end <= window_dw[i])
  cat("gene ",i,"-peaks ",length(sel),"\n")
  cor_p <- vector() # p-value
  cor_r <- vector() # rho
  dist <- vector() # distance to TSS
  gene <- vector() # gene ID
  peak.id <- vector() # peak ID
  if(length(sel)!=0){
    for(j in 1:length(sel)){
      cc <- cor.test(rna[i,],atac[sel[j],],method="spearman")
      cor_p <- c(cor_p,cc$p.value) 
      cor_r <- c(cor_r,cc$estimate)
      if(gene_anno$strand[i]=="-1"){dist <- c(dist,as.numeric(atac_end[sel[j]])-gene_anno$TSS[i])} 
      if(gene_anno$strand[i]=="1"){dist <- c(dist,gene_anno$TSS[i]-as.numeric(atac_start[sel[j]]))}
      gene <- c(gene, gene_anno$ensembl_gene_id[i])
      peak.id <- c(peak.id, rownames(anno_enhancer[sel[j],]))
    }
    cor_output <- data.frame(peak=sel,spearman_pval=cor_p,spearman_rho=cor_r, distance_to_TSS=dist, gene_ID=gene, peak_ID=peak.id)
  }
  if(length(sel)==0){
    cor_output <- NA
  }
  peak_gene_association[[i]] <- cor_output
}

head(peak_gene_association)

#determine the number of contrasts
gene_anno$n_contrasts <- 0

sel <- which((!is.na(peak_gene_association)))

for(i in 1:length(sel)){
  sele <- sel[i]
  gene_anno[sele,]$n_contrasts <- nrow(peak_gene_association[[sele]])
}
head(gene_anno)
a <- gene_anno$n_contrasts
summary(a) #916,242
sum(a)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   28.00   37.00   37.77   47.00  164.00 


#Histogram with number of peaks associated to each gene
require(ggplot2)

axis_size <- element_text(size=14)
data_hist <- data.frame(value=a)

p <- ggplot(data=data_hist, aes(x=value)) +
  geom_histogram(binwidth=1, fill="gray") +
  theme_classic() +
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  labs(title="OCR-Gene associations", 
       x="Number of OCRs", y = "Number of genes")

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/GRN/General/enhancers")
ggsave(plot = p, filename = "hist_gene_ocr_associations.pdf")


##Let's see the signifcance of the enhancer and mrna
adj_peak_gene_association <- lapply(peak_gene_association,FUN=function(X){
  if(is.data.frame(X)){
    X$spearman_adj.pval <- p.adjust(X$spearman_pval,method="fdr") 
  }
  return(X)
})

names(adj_peak_gene_association) <- gene_anno$ensembl_gene_id

spearman_sig_associations_gene_ocr <- lapply(adj_peak_gene_association,FUN=function(X){
  if(is.data.frame(X)){
    return(sum(X$spearman_adj.pval<0.05))
  }
  if(is.na(X)){
    return(X)
  }
})


gene_symbol_sig_spearman <- gene_anno[which(spearman_sig_associations_gene_ocr[]>0),]$gene_id
length(gene_symbol_sig_spearman) #8520

#Save results
data <- do.call("rbind", adj_peak_gene_association)

write.table(data, "OCR_gene_association.txt", sep="\t", dec=".", quote=FALSE, row.names=FALSE, col.names = TRUE)
saveRDS(adj_peak_gene_association, file = "OCR_gene_association.RDS")



#Spearman correlation

#Barplot OCRs (Figure 2C)
data_barplot <- data.frame(peaks=factor(c("0","1","2-5","6-10","11-20",">20"), levels=c("0","1","2-5","6-10","11-20",">20")), value=c(sum(spearman_sig_associations_gene_ocr==0, na.rm=TRUE),sum(spearman_sig_associations_gene_ocr==1, na.rm=TRUE),sum(spearman_sig_associations_gene_ocr>1 & spearman_sig_associations_gene_ocr<6, na.rm=TRUE),sum(spearman_sig_associations_gene_ocr>5 & spearman_sig_associations_gene_ocr<11, na.rm=TRUE),sum(spearman_sig_associations_gene_ocr>10 & spearman_sig_associations_gene_ocr<21, na.rm=TRUE),sum(spearman_sig_associations_gene_ocr>20, na.rm=TRUE)))

##re-size
axis_size <- element_text(size=14)
# Basic barplot
p <- ggplot(data=data_barplot, aes(x=peaks, y=value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=value), vjust=-0.3, size=6)+
  theme_classic() +
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  ylim(0,max(data_barplot$value)+400) +
  labs(title="Significant correlation (Spearman)", 
       x="Number of OCRs", y = "Number of genes")

ggsave(plot = p, filename = "barplot_OCR_gene_association_spearman.pdf")

data_sig <- data[data$spearman_adj.pval < 0.05,]
data_barplot <- data.frame(peaks=factor(c("0","1","2-5","6-10","11-20",">20"), levels=c("0","1","2-5","6-10","11-20",">20")), value=c(sum(table(data_sig$gene_ID)==0, na.rm=TRUE),sum(table(data_sig$gene_ID)==1, na.rm=TRUE),sum(table(data_sig$gene_ID)>1 & table(data_sig$gene_ID)<6, na.rm=TRUE),sum(table(data_sig$gene_ID)>5 & table(data_sig$gene_ID)<11, na.rm=TRUE),sum(table(data_sig$gene_ID)>10 & table(data_sig$gene_ID)<21, na.rm=TRUE),sum(table(data_sig$gene_ID)>20, na.rm=TRUE)))

##re-size
axis_size <- element_text(size=14)
# Basic barplot
p <- ggplot(data=data_barplot, aes(x=peaks, y=value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=value), vjust=-0.3, size=6)+
  theme_classic() +
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  ylim(0,max(data_barplot$value)+400) +
  labs(title="Significant correlation (Spearman)", 
       x="Number of OCRs", y = "Number of genes")

ggsave(plot = p, dpi = 100, filename = "barplot_OCR_gene_sig_association_spearman.pdf")


##Get the list of peaks with and without RNA expression association
list_sig_peaks <- unlist(lapply(adj_peak_gene_association,FUN=function(X){
  if(is.data.frame(X)){
    sel_peak <- X$peak_ID[X$spearman_adj.pval<0.01]
    if(length(sel_peak)>0){return(as.character(sel_peak))}
  }
}))

peaks_sig <- length(unique(data_sig$peak_ID)) #20247

##Distribution of peaks per gene 
aa <- table(list_sig_peaks)

summary(as.numeric(aa))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.0     1.0     1.0     1.7     2.0    19.0  

data_barplot <- data.frame(peaks=factor(c("1","2","3","4","5","6-10",">10"), levels=c("1","2","3","4","5","6-10",">10")), value=c(sum(aa==1, na.rm=TRUE),sum(aa==2, na.rm=TRUE),sum(aa==3, na.rm=TRUE),sum(aa==4, na.rm=TRUE),sum(aa==5, na.rm=TRUE),sum(aa>5 & aa<11, na.rm=TRUE),sum(aa>10, na.rm=TRUE)))

##re-size
axis_size <- element_text(size=14)
# Basic barplot
p <- ggplot(data=data_barplot, aes(x=peaks, y=value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=value), vjust=-0.3, size=6)+
  theme_classic() +
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  ylim(0,max(data_barplot$value)+400) +
  labs(title="Significant correlation (Spearman)", 
       y="Number of enhancer OCRs", x = "Number of genes")

ggsave("barplot_OCR_gene_peaks_sig_association_spearman.pdf")