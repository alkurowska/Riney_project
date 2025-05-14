# Load all ATAC DATA

# Load significant genes coordinated with ChIP 
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load chromatin states 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/Chromatin_states")
states <- read.table("chromatin.states.MM.sorted",sep="\t")
rownames(states) <- paste0(states[,1], "_", states[,2], "_", states[,3])
coln <- read.table("colnames.chromatin.states.MM",sep="\t")
colnames(states) <- coln$V1
dim(states)
states_gr <- states[,c(1:3)]
states <- states[,-c(1:3)]

library(GenomicRanges)
library(ggplot2)

# create GRanges object from the states 
states_gr <- makeGRangesFromDataFrame(states_gr, keep.extra.columns=TRUE)

MGUS_HC_up <- rownames(dea_results[dea_results$sig_MGUS_HC==1,])[!rownames(dea_results[dea_results$sig_MGUS_HC==1,])%in%MGUS_HC]
MGUS_HC_down <- rownames(dea_results[dea_results$sig_MGUS_HC==-1,])[!rownames(dea_results[dea_results$sig_MGUS_HC==-1,])%in%MGUS_HC]
SMM_HC_up <- rownames(dea_results[dea_results$sig_SMM_HC==1,])[!rownames(dea_results[dea_results$sig_SMM_HC==1,])%in%SMM_HC]
SMM_HC_down <- rownames(dea_results[dea_results$sig_SMM_HC==-1,])[!rownames(dea_results[dea_results$sig_SMM_HC==-1,])%in%SMM_HC]
MM_HC_up <- rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])[!rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])%in%MM_HC]
MM_HC_down <- rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])[!rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])%in%MM_HC]

# Make GRanges object from the peaks
library(stringr)
MGUS_HC_up <- as.data.frame(str_split_fixed(MGUS_HC_up, "_", 3))
colnames(MGUS_HC_up) <- c("chr", "start", "end")
MGUS_HC_down <- as.data.frame(str_split_fixed(MGUS_HC_down, "_", 3))
colnames(MGUS_HC_down) <- c("chr", "start", "end")
SMM_HC_up <- as.data.frame(str_split_fixed(SMM_HC_up, "_", 3))
colnames(SMM_HC_up) <- c("chr", "start", "end")
SMM_HC_down <- as.data.frame(str_split_fixed(SMM_HC_down, "_", 3))
colnames(SMM_HC_down) <- c("chr", "start", "end")
MM_HC_up <- as.data.frame(str_split_fixed(MM_HC_up, "_", 3))
colnames(MM_HC_up) <- c("chr", "start", "end")
MM_HC_down <- as.data.frame(str_split_fixed(MM_HC_down, "_", 3))
colnames(MM_HC_down) <- c("chr", "start", "end")

# Make GRanges object from the peaks
MGUS_HC_up <- makeGRangesFromDataFrame(MGUS_HC_up, keep.extra.columns=TRUE)
MGUS_HC_down <- makeGRangesFromDataFrame(MGUS_HC_down, keep.extra.columns=TRUE)
SMM_HC_up <- makeGRangesFromDataFrame(SMM_HC_up, keep.extra.columns=TRUE)
SMM_HC_down <- makeGRangesFromDataFrame(SMM_HC_down, keep.extra.columns=TRUE)
MM_HC_up <- makeGRangesFromDataFrame(MM_HC_up, keep.extra.columns=TRUE)
MM_HC_down <- makeGRangesFromDataFrame(MM_HC_down, keep.extra.columns=TRUE)


# FIND OVERLAPS with states peaks
# Find overlaps
#-------------------------------------

overlaps <- findOverlaps(states_gr, MGUS_HC_up, minoverlap=40) # 22700
MGUS_HC_up <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, MGUS_HC_down, minoverlap=40) # 13738
MGUS_HC_down <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, SMM_HC_up, minoverlap=40) # 35527
SMM_HC_up <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, SMM_HC_down, minoverlap=40) # 23210
SMM_HC_down <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, MM_HC_up, minoverlap=40) # 52465
MM_HC_up <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, MM_HC_down, minoverlap=40) # 43837
MM_HC_down <- states_gr[queryHits(overlaps)]


order.chrom.states <- c("E1","E3","E4","E2","E6","E5","E7","E9","E8","E10","E12","E11")
chr.states.names <- c("Act.Promoter", "Weak.Promoter", "Poised.Promoter", "Strong.Enh1", "Strong.Enh2", "Weak.Enh", "Txn.Transition", "Txn.Elongation", "Weak.Txn", "Heterochr.Rep", "PcG.Rep", "Heterochr")
color_states <- c("#FF0100", "#F6CAE4", "#6E1E8C", "#FFDC64", "#FE8200", "#FEFF54", "#028B64", "#02FF00", "#D3FECD", "#787878", "#AAAAAA", "#F0F0F0")
names(color_states) <- chr.states.names

cell_type <- c("MM", "pb-NBC", "t-NBC", "GCBC", "MBC", "t-PC")
color_cell_type <- c("#8D3C35", "#A6C962", "#525C2A", "#CD2D1F", "#2A61A5", "#B1A3C1")

# Anno plot
# Get the states
MGUS_HC_up_states <- states[names(MGUS_HC_up),]
MGUS_HC_down_states <- states[names(MGUS_HC_down),]
SMM_HC_up_states <- states[names(SMM_HC_up),]
SMM_HC_down_states <- states[names(SMM_HC_down),]
MM_HC_up_states <- states[names(MM_HC_up),]
MM_HC_down_states <- states[names(MM_HC_down),]

metadata <- as.data.frame(coln[-c(1:3),])
colnames(metadata) <- c("sample")
metadata$cell_type <- c("MM", "MM", "MM", "MM", "pb-NBC", "pb-NBC", "pb-NBC", "t-NBC", "t-NBC", "t-NBC", "GCBC", "GCBC", "GCBC", "MBC", "MBC", "MBC", "t-PC", "t-PC", "t-PC")


library(dplyr)
library(tidyverse)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/Chromatin_states_not_coordinated")

list_toPlot <- list(MGUS_HC_up_states, MGUS_HC_down_states, SMM_HC_up_states, SMM_HC_down_states, MM_HC_up_states, MM_HC_down_states)
names(list_toPlot) <- c("MGUS_HC_up", "MGUS_HC_down", "SMM_HC_up", "SMM_HC_down", "MM_HC_up", "MM_HC_down")


for (j in 1:length(list_toPlot)) {
    
    # create a data frame to store the dat
    data_toPlot <- as.data.frame(matrix(nrow=length(order.chrom.states), ncol=0))
    rownames(data_toPlot) <- order.chrom.states

    ploting <- list_toPlot[[j]]
    ploting <- as.data.frame(t(ploting))
    ploting$cell_type <- metadata$cell_type
    
    # Convert to Percentage
    for (i in 1:nrow(ploting)) {
        # get data for sample i
        data <- as.data.frame(t(ploting[i,c(1:(ncol(ploting)-1))]))
        colnames(data) <- "chrom_state"
        data$chrom_state <- factor(data$chrom_state, levels = order.chrom.states)
        data <- data %>% 
                group_by(chrom_state) %>% 
                summarise(freq = n(), .groups = "drop") %>% 
                complete(chrom_state = order.chrom.states, fill = list(freq = 0)) %>%
                mutate(freq = freq/sum(freq)*100)
        for (k in 1:length(order.chrom.states)) {
            if (order.chrom.states[k] %in% data$chrom_state) {
                data_toPlot[k,i] <- data$freq[data$chrom_state == order.chrom.states[k]]
            } else {
                data_toPlot[k,i] <- 0
            }
        }
    }

        colnames(data_toPlot) <- rownames(ploting)

        data_toPlot$chrom_states <- rownames(data_toPlot)
        data_toPlot$chrom_states <- factor(data_toPlot$chrom_states, levels = order.chrom.states)
        data_toPlot$chrom_states_name <- NA 
    for(k in 1:nrow(data_toPlot)){
         data_toPlot$chrom_states_name[k] <- chr.states.names[which(order.chrom.states==data_toPlot$chrom_states[k])]
    }

    # change to a long format
    toPlot_final <- data_toPlot[,c(1:(ncol(data_toPlot)-2),ncol(data_toPlot))]
    rownames(toPlot_final) <- toPlot_final$chrom_states_name
    toPlot_final <- toPlot_final %>% 
        pivot_longer(cols = -chrom_states_name, names_to = "sample", 
        values_to = "freq")
  
    toPlot_final$chrom_states_name <- factor(toPlot_final$chrom_states_name, levels = rev(chr.states.names))
    p <- ggplot(toPlot_final, aes(fill=chrom_states_name, x = sample, y = freq)) + 
    geom_bar(position="fill", stat="identity") +
    # change the order of the x axis
    scale_x_discrete(limits=colnames(data_toPlot)[1:19]) +
    scale_fill_manual(values=color_states) +
    # legend title
    labs(fill="Region") +
    ylab("Percentage") +
    xlab("") + 
    ggtitle("") +
    # change x axis labels 45 degress
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(p, filename = paste0("chrom.states_barplot_", names(list_toPlot)[j], ".pdf"), width = 6, height = 4)
}



# Genes to plot for sankey

MGUS_HC_up <- rownames(dea_results[dea_results$sig_MGUS_HC==1,])[!rownames(dea_results[dea_results$sig_MGUS_HC==1,])%in%MGUS_HC]
MGUS_HC_down <- rownames(dea_results[dea_results$sig_MGUS_HC==-1,])[!rownames(dea_results[dea_results$sig_MGUS_HC==-1,])%in%MGUS_HC]
SMM_HC_up <- rownames(dea_results[dea_results$sig_SMM_HC==1,])[!rownames(dea_results[dea_results$sig_SMM_HC==1,])%in%SMM_HC]
SMM_HC_down <- rownames(dea_results[dea_results$sig_SMM_HC==-1,])[!rownames(dea_results[dea_results$sig_SMM_HC==-1,])%in%SMM_HC]
MM_HC_up <- rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])[!rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])%in%MM_HC]
MM_HC_down <- rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])[!rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])%in%MM_HC]

genes_toPlot <- unique(c(MGUS_HC_up, SMM_HC_up, MM_HC_up, MGUS_HC_down, SMM_HC_down, MM_HC_down))

toPlot <- dea_results[genes_toPlot, c("sig_MGUS_HC","sig_SMM_HC")]
toPlot <- cbind(toPlot, dea_res_noscore[genes_toPlot, c("sig_MM_HC")])

#### PREP THE DATA ####
colnames(toPlot) <- c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC")

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Gene <- rownames(toPlot)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial
genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(MGUS_vs_HC, SMM_vs_HC, MM_vs_HC),
  names_to = "Condition",
  values_to = "Category"
)

# Order conditions
genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC"))

# Plot
png("atac_dea_sanky.png", width=1500, height=1000, res=300)
ggplot(genes_long, aes(
  x = Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
dev.off()
