# Get metadata from mouse
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/")

metadata <- read.csv("metadata.csv")
toPlot <- metadata[,c("Stage", "X.PC.BM.")]

# scale from 0 to 1
# remove NA values
toPlot <- toPlot[!is.na(toPlot$X.PC.BM.),]
toPlot$X.PC.BM. <- (toPlot$X.PC.BM. - min(toPlot$X.PC.BM.)) / (max(toPlot$X.PC.BM.) - min(toPlot$X.PC.BM.))

library(ggplot2)

color_stage <- c("#368716", "#368716",  "#D27D46", "#D92F5E", "#8D4E85") 
names(color_stage) <- c("Control", "HC", "MGUS", "SMM", "MM")

png("PC_boxplot.png", width = 600, height = 600)
ggplot(toPlot, aes(x = Stage, y = X.PC.BM., fill = Stage)) +
  geom_boxplot() +
  scale_fill_manual(values = color_stage) +
  labs(title = "Bone marrow %PC infiltration per stage", x = "Stage", y = "PC BM") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# Plot for human

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")

### STEP 1: LOAD METADATA
load("norm_to_voom.RData")

fish_rna <- fish_rna

# Add entropy and GSVA scores
scores <- read.table("Bulk_samples_scores_all.txt", header=TRUE, sep="")
rownames(scores) <- scores$Sample
table(rownames(fish_rna) == rownames(scores))
fish_rna$gsva <- scores[rownames(fish_rna), "GSVA_low"]

# scale from 0 to 1
toPlot <- fish_rna[,c("Stage", "gsva")]
toPlot$gsva <- (toPlot$gsva - min(toPlot$gsva)) / (max(toPlot$gsva) - min(toPlot$gsva))
toPlot$Stage <- factor(toPlot$Stage, levels = c( "HC", "MGUS", "SMM", "MM"))

png("GSVA_boxplot.png", width = 600, height = 600)
ggplot(toPlot, aes(x = Stage, y = gsva, fill = Stage)) +
  geom_boxplot() +
  scale_fill_manual(values = color_stage) +
  labs(title = "Bone marrow %HC score", x = "Stage", y = "HC proxy") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# without scaling 
toPlot <- fish_rna[,c("Stage", "gsva")]
toPlot$Stage <- factor(toPlot$Stage, levels = c( "HC", "MGUS", "SMM", "MM"))

png("GSVA_boxplot_not_scaled.png", width = 600, height = 600)
ggplot(toPlot, aes(x = Stage, y = gsva, fill = Stage)) +
  geom_boxplot() +
  scale_fill_manual(values = color_stage) +
  labs(title = "Bone marrow %HC score", x = "Stage", y = "HC proxy") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

pdf("GSVA_boxplot_not_scaled_violin.pdf", width = 8, height = 8)
ggplot(toPlot, aes(x = Stage, y = gsva, fill = Stage)) +
  geom_violin() +
  scale_fill_manual(values = color_stage) +
  labs(title = "Bone marrow %HC score", x = "Stage", y = "HC proxy") +
  # size of labels
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))
dev.off()