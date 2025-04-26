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

# Get differntial results of ATAC-Seq
# LOAD DEGs

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
logFC <- as.data.frame(cbind(dea_results$logFC_MGUS_HC, dea_results$logFC_SMM_HC, dea_res_noscore$logFC_MM_HC))
rownames(logFC) <- rownames(dea_results)
colnames(logFC) <- c("MGUS_HC", "SMM_HC", "MM_HC")

# Keep only the peaks that are significant 
logFC <- logFC[all_peaks,]

###################################
# PREPARE PATHWAYS with GENE-SETS #
###################################


# load pathways
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
custom_pathways <- gmtPathways("custom_pathways.gmt")


# MGUS_HC
odds <- data.frame()
chi <- data.frame()

for(i in 1:length(custom_pathways)){
  
  TF <- custom_pathways[[i]]
  
  peaks <- logFC
  # assign the direction
  peaks$direction <- rep("closed",nrow(peaks)) 
  peaks$direction[peaks$MGUS_HC>0] <- "open"
  
  # TFBS
  peaks$TF_binded <- rep("no",nrow(peaks)) 
  peaks[TF,]$TF_binded <- "yes"
  
  x1 <- sum(peaks$direction=="open"&peaks$TF_binded=="yes")
  x2 <- sum(peaks$direction=="closed"&peaks$TF_binded=="yes")
  x3 <- sum(peaks$direction=="open"&peaks$TF_binded=="no")
  x4 <- sum(peaks$direction=="closed"&peaks$TF_binded=="no")
  
  x <- cbind(x1,x3)
  xx <- cbind(x2,x4)

  test.table <- rbind(x, xx)

  colnames(test.table) <- c("TFBS", "no_TFBS")
  rownames(test.table) <- c("open", "closed")

  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(test.table)

  # Extract the odds ratio
  odds_ratio <- fisher_test_result$estimate
  # Extract the p.value
  p.val <- fisher_test_result$p.value

  odds_results <- cbind(odds_ratio, p.val)
  odds <- rbind(odds, odds_results)

  # Perform Chi-Square Test
  chi_square_test_result <- chisq.test(test.table)

  # Extract the X-squared
  X_squared <- chi_square_test_result$statistic
  # Extract the p.value
  X_p.val <- chi_square_test_result$p.val

  chi_results <- cbind(X_squared, X_p.val)
  chi <- rbind(chi, chi_results)
  
}


rownames(odds) <- names(custom_pathways)
rownames(chi) <- names(custom_pathways)

odds$p.adj <- p.adjust(odds$p.val, method = "fdr")
chi$p.adj <- p.adjust(chi$X_p.val, method = "fdr")

# save the results
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/02_OR")
write.table(odds, "MGUS_HC_odds.txt", sep = '\t', col.names = T, row.names = T)
write.table(chi, "MGUS_HC_chi.txt", sep = '\t', col.names = T, row.names = T)




# SMM_HC
odds <- data.frame()
chi <- data.frame()

for(i in 1:length(custom_pathways)){
  
  TF <- custom_pathways[[i]]
  
  peaks <- logFC
  # assign the direction
  peaks$direction <- rep("closed",nrow(peaks)) 
  peaks$direction[peaks$SMM_HC>0] <- "open"
  
  # TFBS
  peaks$TF_binded <- rep("no",nrow(peaks)) 
  peaks[TF,]$TF_binded <- "yes"
  
  x1 <- sum(peaks$direction=="open"&peaks$TF_binded=="yes")
  x2 <- sum(peaks$direction=="closed"&peaks$TF_binded=="yes")
  x3 <- sum(peaks$direction=="open"&peaks$TF_binded=="no")
  x4 <- sum(peaks$direction=="closed"&peaks$TF_binded=="no")
  
  x <- cbind(x1,x3)
  xx <- cbind(x2,x4)

  test.table <- rbind(x, xx)

  colnames(test.table) <- c("TFBS", "no_TFBS")
  rownames(test.table) <- c("open", "closed")

  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(test.table)

  # Extract the odds ratio
  odds_ratio <- fisher_test_result$estimate
  # Extract the p.value
  p.val <- fisher_test_result$p.value

  odds_results <- cbind(odds_ratio, p.val)
  odds <- rbind(odds, odds_results)

  # Perform Chi-Square Test
  chi_square_test_result <- chisq.test(test.table)

  # Extract the X-squared
  X_squared <- chi_square_test_result$statistic
  # Extract the p.value
  X_p.val <- chi_square_test_result$p.val

  chi_results <- cbind(X_squared, X_p.val)
  chi <- rbind(chi, chi_results)
  
}


rownames(odds) <- names(custom_pathways)
rownames(chi) <- names(custom_pathways)

odds$p.adj <- p.adjust(odds$p.val, method = "fdr")
chi$p.adj <- p.adjust(chi$X_p.val, method = "fdr")

# save the results
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/02_OR")
write.table(odds, "SMM_HC_odds.txt", sep = '\t', col.names = T, row.names = T)
write.table(chi, "SMM_HC_chi.txt", sep = '\t', col.names = T, row.names = T)





# MM_HC
odds <- data.frame()
chi <- data.frame()

for(i in 1:length(custom_pathways)){
  
  TF <- custom_pathways[[i]]
  
  peaks <- logFC
  # assign the direction
  peaks$direction <- rep("closed",nrow(peaks)) 
  peaks$direction[peaks$MM_HC>0] <- "open"
  
  # TFBS
  peaks$TF_binded <- rep("no",nrow(peaks)) 
  peaks[TF,]$TF_binded <- "yes"
  
  x1 <- sum(peaks$direction=="open"&peaks$TF_binded=="yes")
  x2 <- sum(peaks$direction=="closed"&peaks$TF_binded=="yes")
  x3 <- sum(peaks$direction=="open"&peaks$TF_binded=="no")
  x4 <- sum(peaks$direction=="closed"&peaks$TF_binded=="no")
  
  x <- cbind(x1,x3)
  xx <- cbind(x2,x4)

  test.table <- rbind(x, xx)

  colnames(test.table) <- c("TFBS", "no_TFBS")
  rownames(test.table) <- c("open", "closed")

  # Perform Fisher's Exact Test
  fisher_test_result <- fisher.test(test.table)

  # Extract the odds ratio
  odds_ratio <- fisher_test_result$estimate
  # Extract the p.value
  p.val <- fisher_test_result$p.value

  odds_results <- cbind(odds_ratio, p.val)
  odds <- rbind(odds, odds_results)

  # Perform Chi-Square Test
  chi_square_test_result <- chisq.test(test.table)

  # Extract the X-squared
  X_squared <- chi_square_test_result$statistic
  # Extract the p.value
  X_p.val <- chi_square_test_result$p.val

  chi_results <- cbind(X_squared, X_p.val)
  chi <- rbind(chi, chi_results)
  
}


rownames(odds) <- names(custom_pathways)
rownames(chi) <- names(custom_pathways)

odds$p.adj <- p.adjust(odds$p.val, method = "fdr")
chi$p.adj <- p.adjust(chi$X_p.val, method = "fdr")

# save the results
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/02_OR")
write.table(odds, "MM_HC_odds.txt", sep = '\t', col.names = T, row.names = T)
write.table(chi, "MM_HC_chi.txt", sep = '\t', col.names = T, row.names = T)