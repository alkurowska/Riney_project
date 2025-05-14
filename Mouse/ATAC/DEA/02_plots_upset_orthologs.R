###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

library(UpSetR)

# Load orthologs 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
orthologs <- read.table("human_mouse_orthologs.txt", sep = "\t", header = TRUE)

# Load mouse peak anno
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Peaks/")
peak_anno <- readRDS("peak_annotation.RDS")
dim(peak_anno) # 51685

# Add orthologs to peak annotation
peak_anno$ortholog <- orthologs$Human_ensembl_gene_id[match(peak_anno$ENSEMBL, orthologs$Mouse_ensembl_gene_id)]
rownames(peak_anno) <- paste0(peak_anno$seqnames, "_", peak_anno$start, "_", peak_anno$end)

# Remove rows with NA orthologs
peak_anno <- peak_anno[!is.na(peak_anno$ortholog),]
head(peak_anno)

# For each ortholog in the peak annotation, find the corresponding human peak annotation
# Get the human peak annotation
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks/")
peak_anno_human <- readRDS("peak_annotation.RDS")
dim(peak_anno_human) # 142136
# remove rows with NA ENSEMBL
peak_anno_human <- peak_anno_human[!is.na(peak_anno_human$ENSEMBL),] # 115485
peak_anno$peaks_human <- NA

for(i in 1:nrow(peak_anno)){
  human_gene <- peak_anno$ortholog[i]
  peaks <- unique(peak_anno_human[peak_anno_human$ENSEMBL == human_gene,]$name)
  # add peaks to the peak_anno separated by comma
  if(length(peaks) > 0){
    peak_anno$peaks_human[i] <- paste(peaks, collapse = ",")
  } else {
    peak_anno$peaks_human[i] <- NA
  }
}

# Remove rows with NA peaks_human
peak_anno <- peak_anno[!is.na(peak_anno$peaks_human),]
dim(peak_anno) # 31108

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA/")
save(peak_anno, file = "peak_anno_orthologs.RData")

# Create a list of peaks pairs between mouse and human
peak_pairs <- list()
for(i in 1:nrow(peak_anno)){
  mouse_peak <- rownames(peak_anno)[i]
  human_peaks <- unlist(strsplit(peak_anno$peaks_human[i], ","))
  peak_pairs[[mouse_peak]] <- human_peaks
}

# change into a data frame each part of the list
# Create a data frame with mouse peaks and their corresponding human peaks
peak_pairs_df <- data.frame()
for(i in 1:length(peak_pairs)){
  mouse_peak <- names(peak_pairs)[i]
  human_peaks <- peak_pairs[[i]]
  for(j in 1:length(human_peaks)){
    peak_pairs_df <- rbind(peak_pairs_df, data.frame(mouse_peak = mouse_peak, human_peak = human_peaks[j]))
  }
}

head(peak_pairs_df) # 31108 peaks
# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/DEA/")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")
rownames(dea_results)
peak_pairs_df <- merge(peak_pairs_df, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", 
                                                     "sig_Trp53_MGUS", "sig_Trp53_MM", 
                                                     "sig_Mmset_MGUS", "sig_Mmset_MM", 
                                                     "sig_MIc_MGUS", "sig_MIc_MM" )], by.x = "mouse_peak", by.y = "row.names")


setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA/")
human_dea <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")
peak_pairs_df <- merge(peak_pairs_df, human_dea[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "human_peak", by.y = "row.names")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
human_dea_no_score <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")
peak_pairs_df$sig_MM_HC <- NA
for (i in 1:nrow(peak_pairs_df)){
  human_peak <- peak_pairs_df$human_peak[i]
  if(human_peak %in% rownames(human_dea_no_score)){
    peak_pairs_df$sig_MM_HC[i] <- human_dea_no_score[human_peak, "sig_MM_HC"]
  }
}

head(peak_pairs_df) # 31108 peaks
dim(peak_pairs_df) # 31108 peaks





