##### READ DATA #####
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
TF_data <- read.table("TFs.txt", header = TRUE)
head(TF_data)
rownames(TF_data) <- make.unique(TF_data$TFs)

# Load libraries
library("RobustRankAggreg")

# Layers
TF_m <- TF_data[,c(4:ncol(TF_data))]
TF_m <- as.matrix(TF_m)
rownames(TF_m) <- rownames(TF_data)

# List of essential TFs
Pioneer <- as.character( read.table("Pioneers.txt",header = F)[,1])
Crisp <- as.character(read.table("Crispr.txt",header = F)[,1])

#Pioneer[!(Pioneer %in% rownames(TF_m))]
#Crisp[!(Crisp %in% rownames(TF_m))]

##### Preparing Transitions ranking #####
# Remove NAs
TF_m <- TF_m[complete.cases(TF_m),]
colnames(TF_m)
TF_MGUS_rank <- TF_m[,c("OR_MGUS_HC","pos_NES_MGUS_HC", "chromVAR_logFC_MGUS_HC")]
TF_SMM_rank <- TF_m[,c("OR_SMM_HC","pos_NES_SMM_HC", "chromVAR_logFC_SMM_HC")]
TF_MM_rank <- TF_m[,c("OR_MM_HC","pos_NES_MM_HC", "chromVAR_logFC_MM_HC")]

# All of the layers
MGUS_list <- list(rownames(TF_MGUS_rank[order(as.numeric(TF_MGUS_rank[,1]), decreasing = T),]),
                rownames(TF_MGUS_rank[order(as.numeric(TF_MGUS_rank[,2]), decreasing = T),]),
                rownames(TF_MGUS_rank[order(as.numeric(TF_MGUS_rank[,3]), decreasing = T),]))

SMM_list <- list(rownames(TF_SMM_rank[order(as.numeric(TF_SMM_rank[,1]), decreasing = T),]),
                rownames(TF_SMM_rank[order(as.numeric(TF_SMM_rank[,2]), decreasing = T),]),
                rownames(TF_SMM_rank[order(as.numeric(TF_MGUS_rank[,3]), decreasing = T),]))

MM_list <- list(rownames(TF_MM_rank[order(as.numeric(TF_MM_rank[,1]), decreasing = T),]),
                rownames(TF_MM_rank[order(as.numeric(TF_MM_rank[,2]), decreasing = T),]),
                rownames(TF_MM_rank[order(as.numeric(TF_MGUS_rank[,3]), decreasing = T),]))



# MGUS 
Rank_MGUS <- aggregateRanks(glist = MGUS_list, N = nrow(TF_data))
rownames(Rank_MGUS)[1:5]  #"JUND.1" "NFE2"   "JUNB.1" "BACH1"  "JDP2" 

# SMM
Rank_SMM <- aggregateRanks(glist = SMM_list, N = nrow(TF_data))
rownames(Rank_SMM)[1:5] # "JUND.1" "NFE2"   "JUNB.1" "BACH1"  "JDP2" 

# MM
Rank_MM <- aggregateRanks(glist = MM_list, N = nrow(TF_data))
rownames(Rank_MM)[1:5] # "JUND.1" "NFE2"   "JUNB.1" "BACH1"  "JDP2"  

###### SIGNIFICANT #######

# p.value < 0.05
Rank_MGUS_sig  <- Rank_MGUS[Rank_MGUS[,"Score"] < 0.05,] 
Rank_SMM_sig  <- Rank_SMM[Rank_SMM[,"Score"] < 0.05,]
Rank_MM_sig  <- Rank_MM[Rank_MM[,"Score"] < 0.05,] 

# PIONEERS in significant
Rank_MGUS_sig[rownames(Rank_MGUS_sig)%in%Pioneer,] #FOXA1
Rank_SMM_sig[rownames(Rank_SMM_sig)%in%Pioneer,]  #FOXA1
Rank_MM_sig[rownames(Rank_MM_sig)%in%Pioneer,]  #FOXA1

# Crisprs in significant 
Rank_MGUS_sig[rownames(Rank_MGUS_sig)%in%Crisp,] # JUND.1 FOS
Rank_SMM_sig[rownames(Rank_SMM_sig)%in%Crisp,] # JUND.1 FOS
Rank_MM_sig[rownames(Rank_MM_sig)%in%Crisp,] # JUND.1 FOS


# Save the data
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking/strict/up")
write.table(Rank_MGUS, file = "Rank_MGUS.txt", sep = "\t", quote = F)
write.table(Rank_SMM, file = "Rank_SMM.txt", sep = "\t", quote = F)
write.table(Rank_MM, file = "Rank_MM.txt", sep = "\t", quote = F)

