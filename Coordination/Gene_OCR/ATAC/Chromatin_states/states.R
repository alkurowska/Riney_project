# Load all ATAC DATA

# Load significant genes coordinated with ChIP 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/ATAC")
all_up <- readRDS("all_up.RDS")
all_down <- readRDS("all_down.RDS")
down_SMM_MM <- readRDS("down_SMM_MM.RDS")
up_SMM_MM <- readRDS("up_SMM_MM.RDS")
unique_MM_up <- readRDS("unique_MM_up.RDS")
unique_MM_down <- readRDS("unique_MM_down.RDS")

# Load chromatin states 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/Chromatin_states")
states <- read.table("chromatin.states.MM.sorted",sep="\t")
rownames(states) <- paste0(states[,1], "_", states[,2], "_", states[,3])
coln <- read.table("colnames.chromatin.states.MM",sep="\t")
colnames(states) <- coln$V1

states_gr <- states[,c(1:3)]
states <- states[,-c(1:3)]


library(GenomicRanges)
library(ggplot2)

# create GRanges object from the states 
states_gr <- makeGRangesFromDataFrame(states_gr, keep.extra.columns=TRUE)

# Make GRanges object from the peaks
library(stringr)
all_up <- as.data.frame(str_split_fixed(all_up, "_", 3))
colnames(all_up) <- c("chr", "start", "end")
all_down <- as.data.frame(str_split_fixed(all_down, "_", 3))
colnames(all_down) <- c("chr", "start", "end")
down_SMM_MM <- as.data.frame(str_split_fixed(down_SMM_MM, "_", 3))
colnames(down_SMM_MM) <- c("chr", "start", "end")
up_SMM_MM <- as.data.frame(str_split_fixed(up_SMM_MM, "_", 3))
colnames(up_SMM_MM) <- c("chr", "start", "end")
unique_MM_up <- as.data.frame(str_split_fixed(unique_MM_up, "_", 3))
colnames(unique_MM_up) <- c("chr", "start", "end")
unique_MM_down <- as.data.frame(str_split_fixed(unique_MM_down, "_", 3))
colnames(unique_MM_down) <- c("chr", "start", "end")

# Make GRanges object from the peaks
all_up <- makeGRangesFromDataFrame(all_up, keep.extra.columns=TRUE)
all_down <- makeGRangesFromDataFrame(all_down, keep.extra.columns=TRUE)
down_SMM_MM <- makeGRangesFromDataFrame(down_SMM_MM, keep.extra.columns=TRUE)
up_SMM_MM <- makeGRangesFromDataFrame(up_SMM_MM, keep.extra.columns=TRUE)
unique_MM_up <- makeGRangesFromDataFrame(unique_MM_up, keep.extra.columns=TRUE)
unique_MM_down <- makeGRangesFromDataFrame(unique_MM_down, keep.extra.columns=TRUE)


# FIND OVERLAPS with states peaks
# Find overlaps
#-------------------------------------

overlaps <- findOverlaps(states_gr, all_up, minoverlap=40) # 965
all_up <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, all_down, minoverlap=40) # 643
all_down <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, down_SMM_MM, minoverlap=40) # 1503
down_SMM_MM <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, up_SMM_MM, minoverlap=40) # 683
up_SMM_MM <- states_gr[queryHits(overlaps)]

overlaps <- findOverlaps(states_gr, unique_MM_up, minoverlap=40) # 1923
unique_MM_up <- states_gr[queryHits(overlaps)]
overlaps <- findOverlaps(states_gr, unique_MM_down, minoverlap=40) # 4912
unique_MM_down <- states_gr[queryHits(overlaps)]
#-------------------------------------


order.chrom.states <- c("E1","E3","E4","E2","E6","E5","E7","E9","E8","E10","E12","E11")
chr.states.names <- c("Act.Promoter", "Weak.Promoter", "Poised.Promoter", "Strong.Enh1", "Strong.Enh2", "Weak.Enh", "Txn.Transition", "Txn.Elongation", "Weak.Txn", "Heterochr.Rep", "PcG.Rep", "Heterochr")
color_states <- c("#FF0100", "#F6CAE4", "#6E1E8C", "#FFDC64", "#FE8200", "#FEFF54", "#028B64", "#02FF00", "#D3FECD", "#787878", "#AAAAAA", "#F0F0F0")
names(color_states) <- chr.states.names

cell_type <- c("MM", "pb-NBC", "t-NBC", "GCBC", "MBC", "t-PC")
color_cell_type <- c("#8D3C35", "#A6C962", "#525C2A", "#CD2D1F", "#2A61A5", "#B1A3C1")

# Anno plot
# Get the states
all_up_states <- states[names(all_up),]
all_down_states <- states[names(all_down),]
down_SMM_MM_states <- states[names(down_SMM_MM),]
up_SMM_MM_states <- states[names(up_SMM_MM),]
unique_MM_up_states <- states[names(unique_MM_up),]
unique_MM_down_states <- states[names(unique_MM_down),]

metadata <- as.data.frame(coln[-c(1:3),])
colnames(metadata) <- c("sample")
metadata$cell_type <- c("MM", "MM", "MM", "MM", "pb-NBC", "pb-NBC", "pb-NBC", "t-NBC", "t-NBC", "t-NBC", "GCBC", "GCBC", "GCBC", "MBC", "MBC", "MBC", "t-PC", "t-PC", "t-PC")


library(dplyr)
library(tidyverse)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/ATAC/Chromatin_states")

list_toPlot <- list(all_up_states, all_down_states, down_SMM_MM_states, up_SMM_MM_states, unique_MM_up_states, unique_MM_down_states)
names(list_toPlot) <- c("all_up", "all_down", "down_SMM_MM", "up_SMM_MM", "unique_MM_up", "unique_MM_down")

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






