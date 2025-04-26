
# Setting working directory
set.seed(1234)

#Loading consensus peaks
### STEP 1: LOAD DATA

#consensus peaks
# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Peaks")
anno <- readRDS("peak_annotation.RDS")
head(anno)
dim(anno[anno$annotation2 == "Distal Intergenic" | anno$annotation2 == "Intron",])

rownames(anno) <- paste0(anno$seqnames, "_", anno$start, "_", anno$end)

# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data <- anno

# add a column with the name of the dataset
data$dataset <- rep(c("ChIP-seq"), c(nrow(anno)))
data$dataset <- factor(data$dataset, levels=c("ChIP-seq"))

# library
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# mutate data to calculate the percentage of each annotation
data <- data %>% 
    group_by(dataset, annotation2) %>% 
    summarise(freq = n()) %>% 
    mutate(freq = freq/sum(freq)*100)

# plot %
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Peaks")
pdf("annotation_barplot.pdf", width=3.25, height=5)
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

