################
### Activity ###
################

# Only up peaks in SMM and MM
# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
load("coordinated_peaks.RData")

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
# Global density plot of logFC
toPlot <- data.frame(dea_results$logFC_SMM_HC, dea_res_noscore$logFC_MM_HC)
colnames(toPlot) <- c("SMM_HC", "MM_HC")
rownames(toPlot) <- rownames(dea_results)

toKeep <- unique(c(SMM_HC, MM_HC))
toPlot <- toPlot[toKeep,]

data <- unique(c(rownames(toPlot)[which(toPlot$SMM_HC>0)],
                   rownames(toPlot)[which(toPlot$MM_HC>0)]))

length(data) #7480 ATAC peaks


# GET associated ChIP Peaks
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Correlations")
peaks_ass <- readRDS("peaks_association.RDS")

# keep only associated with up in ATAC for SMM and MM
data_chip <- unique(peaks_ass[peaks_ass$peak_ID%in%data,]$chip_ID) #7279 ChIP peaks



# Plot violin plot for MGUS, SMM, MM of expression of peaks up peaks in SMM and MM, from nromalzied matrix

# regressing-out gsva and frip for plotting 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
log_cpm_noBatch <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)
colnames(log_cpm_noBatch) <- gsub("X", "", colnames(log_cpm_noBatch))

# Load fish data 

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
load("voom_to_dea.RData")
fish_atac <- fish_atac

up_log_cpm <- log_cpm_noBatch[data,]

library(reshape)
library(dplyr)
library(ggplot2)

up_log_cpm$peak <- rownames(up_log_cpm)
up_log_cpm <- melt(up_log_cpm, id.vars="peak")
up_log_cpm$Stage <- NA
up_log_cpm$Stage <- fish_atac$Stage[match(up_log_cpm$variable, fish_atac$`Sample ATACseq`)]

up_log_cpm$Stage <- factor(up_log_cpm$Stage, levels=c("HC", "MGUS", "SMM", "MM"))

# remove HC samples
up_log_cpm <- up_log_cpm[up_log_cpm$Stage != "HC",]

# Calculate means per stage
mean_values <- up_log_cpm %>%
  group_by(Stage) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  mutate(Stage_numeric = as.numeric(factor(Stage, levels = c("MGUS", "SMM", "MM"))))


setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
pdf("Up_peaks_log_cpm_ATAC.pdf")
# Plotting the violin plot
ggplot(up_log_cpm, aes(x = Stage, y = value, fill = Stage)) +
  geom_violin(trim = FALSE) +
  # Add mean points
  ylim(-2, 10) +
  xlab(NULL) +
  ylab("Log CPM") +
  geom_point(data = mean_values, aes(x = Stage, y = mean_value), shape = 23, size = 2, color = "black", fill = "white") +
  # Connect mean points with a dashed line (no NA issue)
  geom_line(data = mean_values, aes(x = Stage, y = mean_value, group = 1), linetype = "dashed", color = "black") +
  # Colors for stages
  scale_fill_manual(values = c("MGUS" = "#FFBF00", "SMM" = "#FF7F00", "MM" = "#DC143C")) +
  theme_minimal() 
dev.off()


# regressing-out gsva and frip for plotting 
setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix")
log_cpm_noBatch <- read.table("log_cpm_noBatch.txt", header=TRUE, row.names=1)
colnames(log_cpm_noBatch) <- gsub("X", "", colnames(log_cpm_noBatch))

# Load fish data 
load("norm_to_voom.RData")

up_chip_log_cpm <- log_cpm_noBatch[data_chip,]


up_chip_log_cpm$peak <- rownames(up_chip_log_cpm)
up_chip_log_cpm <- melt(up_chip_log_cpm, id.vars="peak")
up_chip_log_cpm$Stage <- NA
up_chip_log_cpm$Stage <- fish_atac$Stage[match(up_chip_log_cpm$variable, fish_atac$`Sample ATACseq`)]


up_chip_log_cpm$Stage <- factor(up_chip_log_cpm$Stage, levels=c("MGUS", "SMM", "MM"))


# Calculate means per stage
mean_values_chip <- up_chip_log_cpm %>%
  group_by(Stage) %>%
  summarise(mean_values_chip = mean(value, na.rm = TRUE)) %>%
  mutate(Stage_numeric = as.numeric(factor(Stage, levels = c("MGUS", "SMM", "MM"))))


setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
pdf("Up_peaks_log_cpm_ChIP.pdf")
# Plotting the violin plot
ggplot(up_chip_log_cpm, aes(x = Stage, y = value, fill = Stage)) +
  geom_violin(trim = FALSE) +
  # Add mean points
  ylim(-2, 10) +
  xlab(NULL) +
  ylab("Log CPM") +
  geom_point(data = mean_values_chip, aes(x = Stage, y = mean_values_chip), shape = 20, size = 2, color = "black", fill = "white") +
  # Connect mean points with a dashed line (no NA issue)
  geom_line(data = mean_values_chip, aes(x = Stage, y = mean_values_chip, group = 1), linetype = "solid", color = "black") +
  # Colors for stages
  scale_fill_manual(values = c("MGUS" = "#FFBF00", "SMM" = "#FF7F00", "MM" = "#DC143C")) +
  theme_minimal() 
dev.off()


# Plot points and lines for ATAC and ChIP in one plot
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")

# Plotting the line/point plot

colnames(mean_values) <- c("Stage", "ATAC", "Stage_numeric")
mean_values$ChIP <- mean_values_chip$mean_values_chip

pdf("Up_peaks_log_cpm_ATAC_ChIP.pdf")
ggplot(mean_values,  aes(x = Stage)) +
  ylim(2, 4) +
  
  # ATAC - Mean Points and Line (Dashed)
  geom_point(aes(y = ATAC, shape = "ATAC"), size = 2, color = "black", fill = "white") +
  geom_line(aes(y = ATAC, linetype = "ATAC", group = 1), color = "black") +
  
  # ChIP - Mean Points and Line (Solid)
  geom_point(aes(y = ChIP, shape = "ChIP"), size = 2, color = "black") +
  geom_line(aes(y = ChIP, linetype = "ChIP", group = 1), color = "black") +
  
  # Labels and Theme
  ylab("Log CPM") +
  xlab(NULL) +
  # Legend customization
  scale_shape_manual(name = NULL, values = c("ATAC" = 23, "ChIP" = 20)) +
  scale_linetype_manual(name = NULL, values = c("ATAC" = "dashed", "ChIP" = "solid")) + 
  theme_minimal()
dev.off()



# PLOT raw
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
load("norm_to_voom.RData")
raw <- as.data.frame(y$counts)
raw <- raw[data,]

library(reshape)
library(dplyr)
library(ggplot2)

raw$peak <- rownames(raw)
raw <- melt(raw, id.vars="peak")
raw$Stage <- NA
raw$Stage <- fish_atac$Stage[match(raw$variable, fish_atac$`Sample ATACseq`)]

raw$Stage <- factor(raw$Stage, levels=c("HC", "MGUS", "SMM", "MM"))

# remove HC samples
raw <- raw[raw$Stage != "HC",]

# Calculate means per stage
mean_values <- raw %>%
  group_by(Stage) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  mutate(Stage_numeric = as.numeric(factor(Stage, levels = c("MGUS", "SMM", "MM"))))

library(ggbreak)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
pdf("Up_peaks_raw_ATAC.pdf")
# Plotting the violin plot
ggplot(raw, aes(x = Stage, y = value, fill = Stage)) +
  geom_violin(trim = FALSE) +
  # Add mean points
  xlab(NULL) +
  ylab("counts") +
  geom_point(data = mean_values, aes(x = Stage, y = mean_value), shape = 23, size = 2, color = "black", fill = "white") +
  # Connect mean points with a dashed line (no NA issue)
  geom_line(data = mean_values, aes(x = Stage, y = mean_value, group = 1), linetype = "dashed", color = "black") +
  # Colors for stages
  scale_fill_manual(values = c("MGUS" = "#FFBF00", "SMM" = "#FF7F00", "MM" = "#DC143C")) +
  theme_minimal() +
  # Use broken y-axis to skip values from 250 to 45000
  # scale_y_continuous(
  #   breaks = c(0, 50, 100, 150, 200, 250, 50000), # Define visible y-axis ticks
  #   labels = c("0", "50", "100", "150", "200", "250", "50000")) +
  scale_y_break(c(250, 48200))  # This creates the visual axis break

dev.off()


setwd("/ibex/user/kurowsaa/Riney_project/ChIP/Count_matrix")
load("norm_to_voom.RData")
raw <- as.data.frame(y$counts)
raw <- raw[data_chip,]

library(reshape)
library(dplyr)
library(ggplot2)

raw$peak <- rownames(raw)
raw <- melt(raw, id.vars="peak")
raw$Stage <- NA
raw$Stage <- fish_chip$Stage[match(raw$variable, fish_chip$`Sample ATACseq`)]

raw$Stage <- factor(raw$Stage, levels=c("MGUS", "SMM", "MM"))

# Calculate means per stage
mean_values_chip <- raw %>%
  group_by(Stage) %>%
  summarise(mean_value = mean(value, na.rm = TRUE)) %>%
  mutate(Stage_numeric = as.numeric(factor(Stage, levels = c("MGUS", "SMM", "MM"))))

library(ggbreak)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")
pdf("Up_peaks_raw_ChIP.pdf")
# Plotting the violin plot
ggplot(raw, aes(x = Stage, y = value, fill = Stage)) +
  geom_violin(trim = FALSE) +
  # Add mean points
  xlab(NULL) +
  ylab("counts") +
  geom_point(data = mean_values_chip, aes(x = Stage, y = mean_value), shape = 20, size = 2, color = "black", fill = "white") +
  # Connect mean points with a dashed line (no NA issue)
  geom_line(data = mean_values_chip, aes(x = Stage, y = mean_value, group = 1), linetype = "solid", color = "black") +
  # Colors for stages
  scale_fill_manual(values = c("MGUS" = "#FFBF00", "SMM" = "#FF7F00", "MM" = "#DC143C")) +
  theme_minimal() +
  # Use broken y-axis to skip values from 250 to 45000
  # scale_y_continuous(
  #   breaks = c(0, 50, 100, 150, 200, 250, 50000), # Define visible y-axis ticks
  #   labels = c("0", "50", "100", "150", "200", "250", "50000")) +
  scale_y_break(c(150, 7344))  # This creates the visual axis break

dev.off()




# Plot points and lines for ATAC and ChIP in one plot
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/ChIP_ATAC/PLOTS")

# Plotting the line/point plot

colnames(mean_values) <- c("Stage", "ATAC", "Stage_numeric")
mean_values$ChIP <- mean_values_chip$mean_value

pdf("Up_peaks_raw_ATAC_ChIP.pdf")
ggplot(mean_values,  aes(x = Stage)) +
  ylim(0, 60) +
  
  # ATAC - Mean Points and Line (Dashed)
  geom_point(aes(y = ATAC, shape = "ATAC"), size = 2, color = "black", fill = "white") +
  geom_line(aes(y = ATAC, linetype = "ATAC", group = 1), color = "black") +
  
  # ChIP - Mean Points and Line (Solid)
  geom_point(aes(y = ChIP, shape = "ChIP"), size = 2, color = "black") +
  geom_line(aes(y = ChIP, linetype = "ChIP", group = 1), color = "black") +
  
  # Labels and Theme
  ylab("counts") +
  xlab(NULL) +
  # Legend customization
  scale_shape_manual(name = NULL, values = c("ATAC" = 23, "ChIP" = 20)) +
  scale_linetype_manual(name = NULL, values = c("ATAC" = "dashed", "ChIP" = "solid")) + 
  theme_minimal()
dev.off()