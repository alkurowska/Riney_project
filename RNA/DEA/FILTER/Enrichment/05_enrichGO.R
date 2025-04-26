#####################################################
###### GENE ONTOLOGY ANALYSIS: CLUSTER PROFILER #####
#####################################################
#Source: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis

set.seed(123)

library(clusterProfiler)
library(org.Hs.eg.db) # Human
#library(org.Mm.eg.db) # Mouse
library(enrichplot)
library(DOSE)
library(ggplot2)



#--------------- 
# LOAD FUNCTIONS
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
source("05_ora.R")
source("05_gse.R")


#--------------- 
# LOAD DATA

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
load("final_genes.RData")

######################
##### MGUS vs HC #####
######################
##### 1. Reading cluster classification ----------
# Create a data frame with gene list and cluster assignment
dea_MGUS_HC <- dea_results[MGUS_HC_filter,]
dea_MGUS_HC$geneID <- rownames(dea_MGUS_HC)
dea_MGUS_HC$cluster <- dea_MGUS_HC$sig_MGUS_HC
dea_MGUS_HC$cluster[dea_MGUS_HC$sig_MGUS_HC == 1] <- "up-regulated"
dea_MGUS_HC$cluster[dea_MGUS_HC$sig_MGUS_HC == -1] <- "down-regulated"

mydata <- cbind(dea_MGUS_HC$geneID, dea_MGUS_HC$cluster)
mydata <- as.data.frame(mydata)
colnames(mydata) <- c("geneID", "cluster")

###### GO terms #######
contrast <- "MGUS_HC"

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
ora(mydata, contrast)


##### 2. Reading enrichment set ----------
# Create a gene vector sorted by logFC
mydata <- dea_results$logFC_MGUS_HC
names(mydata) <- rownames(dea_results)
mydata <- sort(mydata, decreasing = TRUE)

###### GSEA #######
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
gse(mydata, contrast)

# Warning message:
# In fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize,  :
#   For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.


#############################
######### SMM vs HC #########
#############################
##### 1. Reading cluster classification ----------
# Create a data frame with gene list and cluster assignment
dea_SMM_HC <- dea_results[SMM_HC_filter,]
dea_SMM_HC$geneID <- rownames(dea_SMM_HC)
dea_SMM_HC$cluster <- dea_SMM_HC$sig_SMM_HC
dea_SMM_HC$cluster[dea_SMM_HC$sig_SMM_HC == 1] <- "up-regulated"
dea_SMM_HC$cluster[dea_SMM_HC$sig_SMM_HC == -1] <- "down-regulated"

mydata <- cbind(dea_SMM_HC$geneID, dea_SMM_HC$cluster)
mydata <- as.data.frame(mydata)
colnames(mydata) <- c("geneID", "cluster")

###### GO terms #######
contrast <- "SMM_HC"

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
ora(mydata, contrast)


##### 2. Reading enrichment set ----------
# Create a gene vector sorted by logFC
mydata <- dea_results$logFC_SMM_HC
names(mydata) <- rownames(dea_results)

mydata <- sort(mydata, decreasing = TRUE)

###### GSEA #######
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
gse(mydata, contrast)

# Warning message:
# In fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize,  :
#   For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.


#############################
######### MM vs HC #########
#############################


##### 1. Reading cluster classification ----------
# Create a data frame with gene list and cluster assignment
dea_MM_HC <- dea_res_noscore[MM_HC_filter,]
dea_MM_HC$geneID <- rownames(dea_MM_HC)
dea_MM_HC$cluster <- dea_MM_HC$sig_MM_HC
dea_MM_HC$cluster[dea_MM_HC$sig_MM_HC == 1] <- "up-regulated"
dea_MM_HC$cluster[dea_MM_HC$sig_MM_HC == -1] <- "down-regulated"

mydata <- cbind(dea_MM_HC$geneID, dea_MM_HC$cluster)
mydata <- as.data.frame(mydata)
colnames(mydata) <- c("geneID", "cluster")

###### GO terms #######
contrast <- "MM_HC"

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
ora(mydata, contrast)


##### 2. Reading enrichment set ----------
# Create a gene vector sorted by logFC
mydata <- dea_res_noscore$logFC_MM_HC
names(mydata) <- rownames(dea_results)

mydata <- sort(mydata, decreasing = TRUE)

###### GSEA #######
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/Enrichment")
gse(mydata, contrast)

# Warning message:
# In fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize,  :
#   For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.