###########################################
####         ATAC-Seq pipeline         ####
####        TF_motifs_binding.R        ####
####                HUMAN              ####
###########################################


# Setting working directory
set.seed(1234)

# Load libraries
# in bash: mamba activate chromVAR
library(reshape2)
library(scales)
library(fgsea)
library(GenomicRanges)
library(SummarizedExperiment)
library(JASPAR2024)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################
# PREPARE PATHWAYS with GENE-SETS #
###################################


## Getting the required files
# 1- Set of OCR
# Consensus Peaks for the final significant peaks
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_peaks.RData")

final_peaks <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# 2- OCR counts and metadata
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata")
my_counts_matrix <- Counts$counts

peaks <- Counts$annotation
rownames(peaks) <- paste(peaks$Chr, peaks$Start, peaks$End, sep = "_")
rownames(my_counts_matrix) <- rownames(peaks)

# Remove 345_3_S22.sort.rmdup.rmblackls.rmchr.bam
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix) == "345_3_S22.sort.rmdup.rmblackls.rmchr.bam")]

# change colnames of the count data
# Remove everything after "_S"
colnames(my_counts_matrix) <- gsub("_S.*", "", colnames(my_counts_matrix))

# Load metadata

# load metadata
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
fish <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(fish) <- fish[1,]
fish <- fish[-1,]

# Remove FISH data with NA values
fish <- fish[!is.na(fish$`1q amp`),] # 202
dim(fish)

# Remove FISH data with empty RNA seq sampels
fish <- fish[!fish$`Sample RNAseq` == "",] # 197 x 17

my_counts_matrix <- my_counts_matrix[,colnames(my_counts_matrix)%in%fish$`Sample ATACseq`]
dim(my_counts_matrix) # 157817  x  185

fish <- fish[fish$`Sample ATACseq` %in% colnames(my_counts_matrix),]
dim(fish) # 185 x 17

rownames(fish) <- fish$`Sample ATACseq`

# Samples/peaks to keep
my_counts_matrix <- my_counts_matrix[,rownames(fish)]
dim(my_counts_matrix) # 157817  x  185

### Filtered peaks
### STEP 1: LOAD DATA
# Create the GRanges of consensus peaks
peaks <- peaks[final_peaks,]
peak_granges <- makeGRangesFromDataFrame(peaks, keep.extra.columns = T)
peak_granges

my_counts_matrix <- my_counts_matrix[final_peaks,]
table(rownames(my_counts_matrix) == names(peak_granges))



## Creating the object 
table(rownames(fish) == colnames(my_counts_matrix))
fragment_counts <- SummarizedExperiment(assays=list(counts=my_counts_matrix),
                                        rowRanges=peak_granges, colData=fish)

# Getting GC content of peaks

fragment_counts2 <- addGCBias(fragment_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg38)

#Filtering peaks
counts_filtered <- filterPeaks(fragment_counts2, non_overlapping = TRUE)


#Get motifs and what peaks contain motifs
## JASPAR 2024
#https://academic.oup.com/nar/article/46/D1/D252/4616875

#Read the motifs
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
motifs24 <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))

# Save motifs names 

# Save motifs names 
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
names <- as.data.frame(name(motifs24))
save(names, file="JASPAR2024.RData")


# keep only expressed
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
load("norm_to_voom.RData")
genes <- unique(gene_anno[gene_anno$gene_id%in%rownames(d1_norm$counts),]$gene_name)
TFs <- unique(names[,1][names[,1]%in%genes]) #453

motifs24 <- motifs24[which(names[,1]%in%TFs)] #484

# Get motif matches in peaks
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
motif_ix <- matchMotifs(motifs24, counts_filtered, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)

#Save as RData
save(motif_ix, file = "Jaspar_motifs_matches.RData")                        
# Get motif positions within peaks for motifs in peaks 

jaspar_ranges <- matchMotifs(motifs24, peak_granges, genome = "hg38", out = "positions") 
# for(i in 1:length(motifs24)){
#   names(jaspar_ranges)[i] <- unname(name(motifs24[i]))
# }
save(jaspar_ranges, file = "Jaspar_motifs_positions.RData")
##>>>>>>>>>>>>>>>>>>>>> Relate each OCR to motif

jaspar_names <- names(jaspar_ranges)

ocr_jaspar_matrix <- matrix(0,nrow=nrow(peaks),ncol=length(jaspar_names))
colnames(ocr_jaspar_matrix) <- jaspar_names
names(peak_granges) <- paste0(peaks$Chr,"_", peaks$Start, "_", peaks$End)
rownames(ocr_jaspar_matrix) <- names(peak_granges)

for(i in 1:length(jaspar_names)){
  cat(i," - ")
  #select granges for target motif
  jaspar_gr_motif <- jaspar_ranges[[i]]
  #ranges overlaping
  ocr_jaspar_association <- findOverlaps(jaspar_gr_motif,peak_granges,type = "within")
  if(length(ocr_jaspar_association)>0){
    ocr_jaspar_matrix[unique(ocr_jaspar_association@to),i] <- 1    
  }  
}

colnames(ocr_jaspar_matrix) <- jaspar_names

summary(colSums(ocr_jaspar_matrix)==0) #FALSE
write.table(ocr_jaspar_matrix, "ocr_jaspar_matrix.txt", sep = "\t", col.names = T, row.names = T)


# Prepare data for fgsea
peak_TF <- setNames(melt(as.matrix(ocr_jaspar_matrix), varnames = c('peak', 'TF')), c('peak', 'TF', 'values'))
peak_TF <- peak_TF[peak_TF$values != 0,]
dim(peak_TF) # 135,881 pairs of TFs binding to peaks
peak_TF <- data.frame("TermID" = peak_TF$TF, "GeneID" = peak_TF$peak)
head(peak_TF)

# Create a list of TermID with assigned GeneID
custom_pathways <- list()
for(i in 1:length(levels(peak_TF$TermID))){
  TF <- levels(peak_TF$TermID)[i]
  T2G <- peak_TF[peak_TF$TermID == TF,]
  custom_pathways[[TF]] <- c(as.character(unique(T2G$GeneID)))
}
names(custom_pathways) <- levels(peak_TF$TermID)

# Define the output GMT file path
gmt_file <- "custom_pathways.gmt"

# Function to write the pathways to a GMT file
write_gmt <- function(pathways, file) {
  con <- file(file, "wt")
  for (pathway_name in names(pathways)) {
    genes <- pathways[[pathway_name]]
    line <- paste(pathway_name, "Custom pathway", paste(genes, collapse = "\t"), sep = "\t")
    writeLines(line, con)
  }
  close(con)
}

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
# Write the custom pathways to the GMT file
write_gmt(custom_pathways, gmt_file)