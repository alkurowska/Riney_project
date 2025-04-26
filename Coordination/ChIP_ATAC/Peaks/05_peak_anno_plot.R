
# Setting working directory
set.seed(1234)

# library
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

#Loading consensus peaks
### STEP 1: LOAD DATA

#consensus peaks
# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
anno_atac <- readRDS("peak_annotation.RDS")
head(anno_atac)
dim(anno_atac[anno_atac$annotation2 == "Distal Intergenic" | anno_atac$annotation2 == "Intron",])

rownames(anno_atac) <- paste0(anno_atac$seqnames, "_", anno_atac$start, "_", anno_atac$end)

# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data_atac <- anno_atac

# add a column with the name of the dataset
data_atac$dataset <- rep(c("ATAC-seq"), c(nrow(anno_atac)))
data_atac$dataset <- factor(data_atac$dataset, levels=c("ATAC-seq"))

# mutate data to calculate the percentage of each annotation
data_atac <- data_atac %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

#consensus peaks
# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Peaks")
anno_chip <- readRDS("peak_annotation.RDS")
head(anno_chip)
dim(anno_chip[anno_chip$annotation2 == "Distal Intergenic" | anno_chip$annotation2 == "Intron",])

rownames(anno_chip) <- paste0(anno_chip$seqnames, "_", anno_chip$start, "_", anno_chip$end)

# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data_chip <- anno_chip

# add a column with the name of the dataset
data_chip$dataset <- rep(c("ChIP-seq"), c(nrow(anno_chip)))
data_chip$dataset <- factor(data_chip$dataset, levels=c("ChIP-seq"))


# mutate data to calculate the percentage of each annotation
data_chip <- data_chip %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

data <- rbind(data_atac, data_chip)

# plot %
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/Peaks")
pdf("annotation_barplot_all.pdf", width=4, height=5)
ggplot(data, aes(fill=annotation2, y=freq, x=`dataset`)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T) +
    # legend title
    labs(fill="Region") +
    ylab("Percentage") +
    xlab("") + 
    ggtitle("Peaks' annotation") +
    # change x axis labels 45 degress
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Load significant genes 
# GET associated ChIP Peaks
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
peaks_ass <- readRDS("peaks_association.RDS")

# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data_atac <- anno_atac[unique(peaks_ass$peak_ID),]

# add a column with the name of the dataset
data_atac$dataset <- rep(c("ATAC-seq"), c(nrow(anno_atac[unique(peaks_ass$peak_ID),])))
data_atac$dataset <- factor(data_atac$dataset, levels=c("ATAC-seq"))

# mutate data to calculate the percentage of each annotation
data_atac <- data_atac %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)


# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data_chip <- anno_chip[unique(peaks_ass$chip_ID),]

# add a column with the name of the dataset
data_chip$dataset <- rep(c("ChIP-seq"), c(nrow(anno_chip[unique(peaks_ass$chip_ID),])))
data_chip$dataset <- factor(data_chip$dataset, levels=c("ChIP-seq"))


# mutate data to calculate the percentage of each annotation
data_chip <- data_chip %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

data <- rbind(data_atac, data_chip)

# plot %
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/Peaks")
pdf("annotation_barplot_coordinated.pdf", width=4, height=5)
ggplot(data, aes(fill=annotation2, y=freq, x=`dataset`)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T) +
    # legend title
    labs(fill="Region") +
    ylab("Percentage") +
    xlab("") + 
    ggtitle("Peaks' annotation") +
    # change x axis labels 45 degress
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()