
# Setting working directory
set.seed(1234)

# LOAD Differential Results
# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")

#Loading consensus peaks
### STEP 1: LOAD DATA

#consensus peaks
# Include peaks annotation
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Peaks")
anno <- readRDS("peak_annotation.RDS")
head(anno)
dim(anno[anno$annotation2 == "Distal Intergenic" | anno$annotation2 == "Intron",])

rownames(anno) <- paste0(anno$seqnames, "_", anno$start, "_", anno$end)

# HC -> MGUS
MGUS_HC <- anno[MGUS_HC_filter,]

SMM_HC <- anno[SMM_HC_filter,]

MM_HC <- anno[MM_HC_filter,]



# plot a bar plot of the annotation % of the peaks for all of the peaks "anno", and other contrasts 
# create a dataset
data <- rbind(anno, MGUS_HC, SMM_HC, MM_HC)

# add a column with the name of the dataset
data$dataset <- rep(c("ATAC-seq", "MGUS vs HC", "SMM vs HC", "MM vs HC"), c(nrow(anno), nrow(MGUS_HC), nrow(SMM_HC), nrow(MM_HC)))
data$dataset <- factor(data$dataset, levels=c("ATAC-seq", "MGUS vs HC", "SMM vs HC", "MM vs HC"))
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
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/PLOTS")
pdf("annotation_barplot.pdf", width=5, height=5)
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
