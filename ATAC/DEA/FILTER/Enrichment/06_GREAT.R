#####################################################
######     GENE ONTOLOGY ANALYSIS:     GREAT    #####
#####################################################

### STEP 1: LOAD DATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")

##############################
######### MGUS vs HC #########
##############################

##### 1. Reading cluster classification ----------
### STEP 1: LOAD DATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header = T, sep = "\t")
dea_MGUS_HC <- dea_results[MGUS_HC_filter,]
dea_MGUS_HC$geneID <- rownames(dea_MGUS_HC)
dea_MGUS_HC$cluster <- dea_MGUS_HC$sig_MGUS_HC
dea_MGUS_HC$cluster[dea_MGUS_HC$sig_MGUS_HC == 1] <- "up-regulated"
dea_MGUS_HC$cluster[dea_MGUS_HC$sig_MGUS_HC == -1] <- "down-regulated"

mydata <- dea_MGUS_HC[,c("geneID", "cluster")]

# Consensus Peaks 
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# Make GRanges object
library(GenomicRanges)
# Split row names into 3 columns (chr, start, end)
mydata <- cbind(mydata, data.frame(do.call('rbind', strsplit(rownames(mydata), "_"))))
colnames(mydata) <- c("geneID", "cluster", "chr", "start", "end")
great_mydata <- makeGRangesFromDataFrame(mydata, keep.extra.columns = T)

#### 2. The Genomic Regions Enrichment of Annotations Tool -------

library(rGREAT)

cluster_1 <- great_mydata[great_mydata$cluster == "down-regulated",]
cluster_2 <- great_mydata[great_mydata$cluster == "up-regulated",]

# Biological Process 
greatJob_1 <- great(cluster_1, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)
greatJob_2 <- great(cluster_2, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/Enrichment")
save(greatJob_1, greatJob_2, file = "MGUS_HC_greatJob.RData")

great_1 = getEnrichmentTables(greatJob_1)
great_2 = getEnrichmentTables(greatJob_2)

save(great_1, great_2, file = "MGUS_HC_GREAT_ENRICH_clusters.RData")

great_1_BP <- great_1[order(great_1$p_adjust, decreasing = F),]
great_2_BP <- great_2[order(great_2$p_adjust, decreasing = F),]

great_1_BP$cluster <- rep(1, nrow(great_1_BP))
great_2_BP$cluster <- rep(2, nrow(great_2_BP))

great_1_BP <- great_1_BP[great_1_BP$p_adjust <= 0.05,]
great_2_BP <- great_2_BP[great_2_BP$p_adjust <= 0.05,]

great_1_BP_top <- great_1_BP[1:50,]
great_2_BP_top <- great_2_BP[1:50,]

great_1_BP_top <- great_1_BP_top[order(great_1_BP_top$fold_enrichment, decreasing = T),]
great_2_BP_top <- great_2_BP_top[order(great_2_BP_top$fold_enrichment, decreasing = T),]

# Remove NA
great_1_BP_top <- great_1_BP_top[!is.na(great_1_BP_top$cluster),]
great_2_BP_top <- great_2_BP_top[!is.na(great_2_BP_top$cluster),]

#top 10
toPlot_BP <- rbind(great_1_BP_top[1:10,],
                   great_2_BP_top[1:10,])

rownames(toPlot_BP) <- 1:nrow(toPlot_BP)
# change 1 in cluster column into "down-regulated" and 2 into "up-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 1] <- "down-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 2] <- "up-regulated"

library(ggplot2)
BP <- ggplot(toPlot_BP, 
             aes(as.factor(cluster), y=reorder(description,as.factor(cluster)), size = fold_enrichment,
                 colour = p_adjust)) +
  geom_point() +
  ylab(NULL) +
  xlab(NULL) + 
  ggtitle("Biological Processes") +
  theme_bw() +
  theme(legend.title=element_text(size=8)) +
  labs(size = "foldEnrichment",
       color = "padjust")

ggsave(plot = BP, width = 10, height = 8, dpi = 400, filename = "MGUS_HC_GREAT_GO_BP_dotplot.jpg")

write.table(great_1, "MGUS_HC_GREAT_GO_BP_down_regulated.txt", sep = "\t", quote = F)
write.table(great_2, "MGUS_HC_GREAT_GO_BP_up_regulated.txt", sep = "\t", quote = F)


#############################
######### SMM vs HC #########
#############################

##### 1. Reading cluster classification ----------
### STEP 1: LOAD DATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header = T, sep = "\t")
dea_SMM_HC <- dea_results[SMM_HC_filter,]
dea_SMM_HC$geneID <- rownames(dea_SMM_HC)
dea_SMM_HC$cluster <- dea_SMM_HC$sig_SMM_HC
dea_SMM_HC$cluster[dea_SMM_HC$sig_SMM_HC == 1] <- "up-regulated"
dea_SMM_HC$cluster[dea_SMM_HC$sig_SMM_HC == -1] <- "down-regulated"

mydata <- dea_SMM_HC[,c("geneID", "cluster")]

# Consensus Peaks 
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# Make GRanges object
library(GenomicRanges)
# Split row names into 3 columns (chr, start, end)
mydata <- cbind(mydata, data.frame(do.call('rbind', strsplit(rownames(mydata), "_"))))
colnames(mydata) <- c("geneID", "cluster", "chr", "start", "end")
great_mydata <- makeGRangesFromDataFrame(mydata, keep.extra.columns = T)

#### 2. The Genomic Regions Enrichment of Annotations Tool -------

library(rGREAT)

cluster_1 <- great_mydata[great_mydata$cluster == "down-regulated",]
cluster_2 <- great_mydata[great_mydata$cluster == "up-regulated",]

# Biological Process 
greatJob_1 <- great(cluster_1, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)
greatJob_2 <- great(cluster_2, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/Enrichment")
save(greatJob_1, greatJob_2, file = "SMM_HC_greatJob.RData")

great_1 = getEnrichmentTables(greatJob_1)
great_2 = getEnrichmentTables(greatJob_2)

save(great_1, great_2, file = "SMM_HC_GREAT_ENRICH_clusters.RData")

great_1_BP <- great_1[order(great_1$p_adjust, decreasing = F),]
great_2_BP <- great_2[order(great_2$p_adjust, decreasing = F),]

great_1_BP$cluster <- rep(1, nrow(great_1_BP))
great_2_BP$cluster <- rep(2, nrow(great_2_BP))

great_1_BP <- great_1_BP[great_1_BP$p_adjust <= 0.05,]
great_2_BP <- great_2_BP[great_2_BP$p_adjust <= 0.05,]

great_1_BP_top <- great_1_BP[1:50,]
great_2_BP_top <- great_2_BP[1:50,]

great_1_BP_top <- great_1_BP_top[order(great_1_BP_top$fold_enrichment, decreasing = T),]
great_2_BP_top <- great_2_BP_top[order(great_2_BP_top$fold_enrichment, decreasing = T),]

# Remove NA
great_1_BP_top <- great_1_BP_top[!is.na(great_1_BP_top$cluster),]
great_2_BP_top <- great_2_BP_top[!is.na(great_2_BP_top$cluster),]

#top 10
toPlot_BP <- rbind(great_1_BP_top[1:10,],
                   great_2_BP_top[1:10,])

rownames(toPlot_BP) <- 1:nrow(toPlot_BP)
# change 1 in cluster column into "down-regulated" and 2 into "up-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 1] <- "down-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 2] <- "up-regulated"

library(ggplot2)
BP <- ggplot(toPlot_BP, 
             aes(as.factor(cluster), y=reorder(description,as.factor(cluster)), size = fold_enrichment,
                 colour = p_adjust)) +
  geom_point() +
  ylab(NULL) +
  xlab(NULL) + 
  ggtitle("Biological Processes") +
  theme_bw() +
  theme(legend.title=element_text(size=8)) +
  labs(size = "foldEnrichment",
       color = "padjust")

ggsave(plot = BP, width = 10, height = 8, dpi = 400, filename = "SMM_HC_GREAT_GO_BP_dotplot.jpg")

write.table(great_1, "SMM_HC_GREAT_GO_BP_down_regulated.txt", sep = "\t", quote = F)
write.table(great_2, "SMM_HC_GREAT_GO_BP_up_regulated.txt", sep = "\t", quote = F)




#############################
######### SMM vs HC #########
#############################

##### 1. Reading cluster classification ----------
### STEP 1: LOAD DATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header = T, sep = "\t")
dea_MM_HC <- dea_res_noscore[MM_HC_filter,]
dea_MM_HC$geneID <- rownames(dea_MM_HC)
dea_MM_HC$cluster <- dea_MM_HC$sig_MM_HC
dea_MM_HC$cluster[dea_MM_HC$sig_MM_HC == 1] <- "up-regulated"
dea_MM_HC$cluster[dea_MM_HC$sig_MM_HC == -1] <- "down-regulated"

mydata <- dea_SMM_HC[,c("geneID", "cluster")]

# Consensus Peaks 
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# Make GRanges object
library(GenomicRanges)
# Split row names into 3 columns (chr, start, end)
mydata <- cbind(mydata, data.frame(do.call('rbind', strsplit(rownames(mydata), "_"))))
colnames(mydata) <- c("geneID", "cluster", "chr", "start", "end")
great_mydata <- makeGRangesFromDataFrame(mydata, keep.extra.columns = T)

#### 2. The Genomic Regions Enrichment of Annotations Tool -------

library(rGREAT)

cluster_1 <- great_mydata[great_mydata$cluster == "down-regulated",]
cluster_2 <- great_mydata[great_mydata$cluster == "up-regulated",]

# Biological Process 
greatJob_1 <- great(cluster_1, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)
greatJob_2 <- great(cluster_2, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene", basal_upstream = 1.0, basal_downstream = 1.0)

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/Enrichment")
save(greatJob_1, greatJob_2, file = "MM_HC_greatJob.RData")

great_1 = getEnrichmentTables(greatJob_1)
great_2 = getEnrichmentTables(greatJob_2)

save(great_1, great_2, file = "MM_HC_GREAT_ENRICH_clusters.RData")

great_1_BP <- great_1[order(great_1$p_adjust, decreasing = F),]
great_2_BP <- great_2[order(great_2$p_adjust, decreasing = F),]

great_1_BP$cluster <- rep(1, nrow(great_1_BP))
great_2_BP$cluster <- rep(2, nrow(great_2_BP))

great_1_BP <- great_1_BP[great_1_BP$p_adjust <= 0.05,]
great_2_BP <- great_2_BP[great_2_BP$p_adjust <= 0.05,]

great_1_BP_top <- great_1_BP[1:50,]
great_2_BP_top <- great_2_BP[1:50,]

great_1_BP_top <- great_1_BP_top[order(great_1_BP_top$fold_enrichment, decreasing = T),]
great_2_BP_top <- great_2_BP_top[order(great_2_BP_top$fold_enrichment, decreasing = T),]

# Remove NA
great_1_BP_top <- great_1_BP_top[!is.na(great_1_BP_top$cluster),]
great_2_BP_top <- great_2_BP_top[!is.na(great_2_BP_top$cluster),]

#top 10
toPlot_BP <- rbind(great_1_BP_top[1:10,],
                   great_2_BP_top[1:10,])

rownames(toPlot_BP) <- 1:nrow(toPlot_BP)
# change 1 in cluster column into "down-regulated" and 2 into "up-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 1] <- "down-regulated"
toPlot_BP$cluster[toPlot_BP$cluster == 2] <- "up-regulated"

library(ggplot2)
BP <- ggplot(toPlot_BP, 
             aes(as.factor(cluster), y=reorder(description,as.factor(cluster)), size = fold_enrichment,
                 colour = p_adjust)) +
  geom_point() +
  ylab(NULL) +
  xlab(NULL) + 
  ggtitle("Biological Processes") +
  theme_bw() +
  theme(legend.title=element_text(size=8)) +
  labs(size = "foldEnrichment",
       color = "padjust")

ggsave(plot = BP, width = 10, height = 8, dpi = 400, filename = "MM_HC_GREAT_GO_BP_dotplot.jpg")

write.table(great_1, "MM_HC_GREAT_GO_BP_down_regulated.txt", sep = "\t", quote = F)
write.table(great_2, "MM_HC_GREAT_GO_BP_up_regulated.txt", sep = "\t", quote = F)
