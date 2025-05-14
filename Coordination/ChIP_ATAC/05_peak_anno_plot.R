
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

#finalpeaks
# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")

anno_MGUS <- anno_atac[MGUS_HC,]
anno_MGUS$dataset <- rep(c("MGUS_HC"), c(nrow(anno_MGUS)))
anno_MGUS$dataset <- factor(anno_MGUS$dataset, levels=c("MGUS_HC"))
anno_MGUS <- anno_MGUS %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

anno_SMM <- anno_atac[SMM_HC,]
anno_SMM$dataset <- rep(c("SMM_HC"), c(nrow(anno_SMM)))
anno_SMM$dataset <- factor(anno_SMM$dataset, levels=c("SMM_HC"))
anno_SMM <- anno_SMM %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

anno_MM <- anno_atac[MM_HC,]
anno_MM$dataset <- rep(c("MM_HC"), c(nrow(anno_MM)))
anno_MM$dataset <- factor(anno_MM$dataset, levels=c("MM_HC"))
anno_MM <- anno_MM %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

data <- rbind(data_atac, anno_MGUS, anno_SMM, anno_MM)

# plot %

pdf("annotation_barplot_differential.pdf", width=5, height=4)
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
