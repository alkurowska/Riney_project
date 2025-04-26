###########################################
####         ATAC-Seq pipeline         ####
####        GSEA_TF_chromatin.R        ####
####                HUMAN              ####
###########################################


# Setting working directory
set.seed(1234)

# Load libraries
library(reshape2)
library(scales)
library(fgsea)


############################
# PREPARE PRE-RANKED PEAKS #
############################


# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_peaks.RData")

all_peaks <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# create a named vector of rank peaks by logFC values
logFC <- dea_res_noscore$logFC_MM_HC
names(logFC) <- rownames(dea_res_noscore)

ranked_logFC <- sort(logFC, decreasing = T)
head(ranked_logFC)

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/all_peaks/00_custom_pathways")
# load pathways
custom_pathways <- gmtPathways("custom_pathways.gmt")


library(doParallel)
# Register the parallel backend
numCores <- detectCores() - 1  # Use one less than the number of available cores
cl <- makeCluster(numCores)
registerDoParallel(cl)

# absolute
res <- foreach(i = 1:length(custom_pathways), .packages = c("fgsea"), .combine = rbind) %dopar% {
    # Gene set enrichment of TFs in peaks with a personalized gene set
    fgseaRes <- fgsea::fgsea(pathways = custom_pathways[i], scoreType = "pos",
                           stats = ranked_logFC[ranked_logFC>0], minSize = 10, maxSize = length(ranked_logFC) - 1,
                           nPermSimple = 100000)
}

# p.adj
res$p.adjust <- p.adjust(res$pval, method = "fdr")

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/all_peaks/01_ES/MM_HC")
save(res, file = "MM_HC_GSEA_TF_pos.RData")

# Stop the cluster
stopCluster(cl)
registerDoSEQ()  # Return to sequential processing

load("MM_HC_GSEA_TF_pos.RData")
res <- as.data.frame(res)
res <- res[,c("pathway", "ES", "NES", "pval", "p.adjust")]
write.table(res, "MM_HC_GSEA_TF_pos.txt", sep = "\t", col.names = T, row.names = F)
