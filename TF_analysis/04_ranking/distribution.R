# correlation
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
TF_data <- read.table("TFs.txt", header = TRUE)
head(TF_data)

library(ggrepel)
library(ggplot2)

colnames(TF_data)

# Plot with NES pos and NES neg
toPlot <- data.frame(TF_data$pos_NES_MGUS_HC, TF_data$pos_NES_SMM_HC, TF_data$pos_NES_MM_HC,
                    TF_data$neg_NES_MGUS_HC, TF_data$neg_NES_SMM_HC, TF_data$neg_NES_MM_HC,
                                        TF_data$abs_NES_MGUS_HC, TF_data$abs_NES_SMM_HC, TF_data$abs_NES_MM_HC)
colnames(toPlot) <- c("MGUS_HC_NES_pos", "SMM_HC_NES_pos", "MM_HC_NES_pos",
                      "MGUS_HC_NES_neg", "SMM_HC_NES_neg", "MM_HC_NES_neg",
                      "MGUS_HC_NES_abs", "SMM_HC_NES_abs", "MM_HC_NES_abs")

# Add rownames, if duplicated add ".number"
rownames(toPlot) <- make.unique(TF_data$TFs)


Transition <- sapply(strsplit(colnames(toPlot),"_"), `[`, 1)
Test <- gsub("MGUS_HC_", "", colnames(toPlot)) 
Test <- gsub("SMM_HC_", "", Test)
Test <- gsub("MM_HC_", "", Test)

metadata <- data.frame(Transition, Test)
metadata$Test <- factor(metadata$Test, levels = c(
                        "NES_pos",
                        "NES_neg",
                        "NES_abs"))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

color_trans <- c("#D27D46", "#D92F5E", "#8D4E85") 
names(color_trans) <- c("MGUS", "SMM", "MM")


color_test <- c( "darkslategray1", "darkturquoise", 
                "deepskyblue3")
names(color_test) <- unique(Test)


library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Function for density plots
density_plot_NES <- function(toPlot, stage_name) {
    p <- ggplot2::ggplot(toPlot, aes(x = NES_abs, fill = "NES_abs")) +
    geom_density(alpha = 0.9) +
    geom_density(aes(x = NES_pos, fill = "NES_pos"), alpha = 0.9) +
    geom_density(aes(x = NES_neg, fill = "NES_neg"), alpha = 0.9) +
    xlab("NES") + 
    ylab("Density") +
    # dash line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("NES_abs" = "deepskyblue3", 
                                "NES_pos" = "darkslategray1", 
                                "NES_neg" = "darkturquoise")) +
    guides(fill = guide_legend(title = "Score type")) +
    labs(color = "") +
    ggtitle(paste0("Density of NES scores - ", stage_name)) +
    theme_minimal()
    
    ggsave(p, filename = paste0(stage_name, "_density_NES.png"), width = 10, height = 10, dpi = 300)
    return(p)
}

# Create individual plots for each stage
contrasts <- c("MGUS_HC", "SMM_HC", "MM_HC")
plots_list <- list()

for(contrast in contrasts) {
    stage_data <- data.frame(
        NES_abs = toPlot[,paste0(contrast, "_NES_abs")],
        NES_pos = toPlot[,paste0(contrast, "_NES_pos")],
        NES_neg = toPlot[,paste0(contrast, "_NES_neg")]
    )
    plots_list[[contrast]] <- density_plot_NES(stage_data, contrast)
}

# Combined plot
p_combined <- ggplot() +
    # MGUS
    geom_density(data = data.frame(score = toPlot$MGUS_HC_NES_abs), aes(x = score, fill = "MGUS_abs"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$MGUS_HC_NES_pos), aes(x = score, fill = "MGUS_pos"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$MGUS_HC_NES_neg), aes(x = score, fill = "MGUS_neg"), alpha = 0.3) +
    # SMM
    geom_density(data = data.frame(score = toPlot$SMM_HC_NES_abs), aes(x = score, fill = "SMM_abs"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$SMM_HC_NES_pos), aes(x = score, fill = "SMM_pos"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$SMM_HC_NES_neg), aes(x = score, fill = "SMM_neg"), alpha = 0.3) +
    # MM
    geom_density(data = data.frame(score = toPlot$MM_HC_NES_abs), aes(x = score, fill = "MM_abs"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$MM_HC_NES_pos), aes(x = score, fill = "MM_pos"), alpha = 0.3) +
    geom_density(data = data.frame(score = toPlot$MM_HC_NES_neg), aes(x = score, fill = "MM_neg"), alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(
        "MGUS_abs" = "#E0A47C", "MGUS_pos" = "#D27D46", "MGUS_neg" = "#A65728",
        "SMM_abs" = "#E56B8E", "SMM_pos" = "#D92F5E", "SMM_neg" = "#AD2548",
        "MM_abs" = "#B07DAB", "MM_pos" = "#8D4E85", "MM_neg" = "#713E6A"
    )) +
    theme_minimal() +
    labs(x = "NES", y = "Density", title = "Combined NES score distributions") +
    facet_wrap(~"Stage", scales = "free")

ggsave("combined_NES_density.png", p_combined, width = 15, height = 8, dpi = 300)


# Calculate summary statistics - corrected version
summary_stats <- data.frame()

for(stage in c("MGUS", "SMM", "MM")) {
    # For NES_abs
    abs_stats <- data.frame(
        Stage = stage,
        Score_type = "NES_abs",
        Mean = mean(toPlot[,paste0(stage, "_HC_NES_abs")], na.rm = TRUE),
        Median = median(toPlot[,paste0(stage, "_HC_NES_abs")], na.rm = TRUE),
        SD = sd(toPlot[,paste0(stage, "_HC_NES_abs")], na.rm = TRUE)
    )
    
    # For NES_pos
    pos_stats <- data.frame(
        Stage = stage,
        Score_type = "NES_pos",
        Mean = mean(toPlot[,paste0(stage, "_HC_NES_pos")], na.rm = TRUE),
        Median = median(toPlot[,paste0(stage, "_HC_NES_pos")], na.rm = TRUE),
        SD = sd(toPlot[,paste0(stage, "_HC_NES_pos")], na.rm = TRUE)
    )
    
    # For NES_neg
    neg_stats <- data.frame(
        Stage = stage,
        Score_type = "NES_neg",
        Mean = mean(toPlot[,paste0(stage, "_HC_NES_neg")], na.rm = TRUE),
        Median = median(toPlot[,paste0(stage, "_HC_NES_neg")], na.rm = TRUE),
        SD = sd(toPlot[,paste0(stage, "_HC_NES_neg")], na.rm = TRUE)
    )
    
    summary_stats <- rbind(summary_stats, abs_stats, pos_stats, neg_stats)
}

# Add number of TFs for each category
summary_stats$N_TFs <- NA
for(i in 1:nrow(summary_stats)) {
    stage <- summary_stats$Stage[i]
    score_type <- summary_stats$Score_type[i]
    col_name <- paste0(stage, "_HC_", score_type)
    summary_stats$N_TFs[i] <- sum(!is.na(toPlot[,col_name]))
}

write.table(summary_stats, "NES_summary_stats.txt", sep="\t", quote=FALSE, row.names=FALSE)




library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Check for duplicated TFs
duplicated_tfs <- TF_data$TFs[duplicated(TF_data$TFs)]
print("Duplicated TFs:")
print(duplicated_tfs)

# Make unique names for duplicated TFs
TF_data$TF_unique <- make.unique(TF_data$TFs)

# Create a mapping of original to unique names for reference
tf_mapping <- data.frame(
    Original = TF_data$TFs,
    Unique = TF_data$TF_unique
)
write.table(tf_mapping, "TF_name_mapping.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Define color schemes
color_test <- c("darkslategray1", "darkturquoise", "deepskyblue3")
names(color_test) <- c("NES_pos", "NES_neg", "NES_abs")

# Function to create plot for one stage
plot_stage_NES <- function(stage) {
    # Get data for this stage
    data <- data.frame(
        TF = TF_data$TF_unique,  # Use unique names
        Original_TF = TF_data$TFs,  # Keep original names for reference
        NES_pos = TF_data[[paste0("pos_NES_", stage)]],
        NES_neg = TF_data[[paste0("neg_NES_", stage)]],
        NES_abs = TF_data[[paste0("abs_NES_", stage)]]
    )
    
    # Convert to long format
    data_long <- pivot_longer(data, 
                            cols = c(NES_pos, NES_neg, NES_abs),
                            names_to = "Score_type",
                            values_to = "NES")
    
    # Order TFs by absolute NES score
    tf_order <- data %>%
        arrange(desc(abs(NES_abs))) %>%
        pull(TF)
    
    data_long$TF <- factor(data_long$TF, levels = tf_order)
    
    # Identify duplicated TFs for highlighting
    duplicated_tfs <- unique(data$Original_TF[duplicated(data$Original_TF)])
    is_duplicated <- data_long$TF %in% data$TF[data$Original_TF %in% duplicated_tfs]
    
    # Create plot
    p <- ggplot(data_long, aes(x = TF, y = NES, fill = Score_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = color_test) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       color = ifelse(levels(data_long$TF) %in% 
                                                    data$TF[data$Original_TF %in% duplicated_tfs],
                                                    "red", "black")),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        labs(x = "Transcription Factors", 
             y = "NES Score",
             title = paste("NES Scores Comparison -", stage),
             subtitle = "Duplicated TFs shown in red") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    
    # Save plot
    ggsave(paste0("NES_comparison_", stage, ".png"), 
           p, width = 20, height = 8, dpi = 300)
    
    return(p)
}

# Create plots for each stage
stages <- c("MGUS_HC", "SMM_HC", "MM_HC")
plots <- lapply(stages, plot_stage_NES)

# Create a table with all scores
all_scores <- data.frame(
    TF = TF_data$TF_unique,
    Original_TF = TF_data$TFs
)

for(stage in stages) {
    all_scores[[paste0(stage, "_pos")]] <- TF_data[[paste0("pos_NES_", stage)]]
    all_scores[[paste0(stage, "_neg")]] <- TF_data[[paste0("neg_NES_", stage)]]
    all_scores[[paste0(stage, "_abs")]] <- TF_data[[paste0("abs_NES_", stage)]]
}

# Order by absolute values in MM_HC
all_scores <- all_scores[order(abs(all_scores$MM_HC_abs), decreasing = TRUE),]

# Write table with both original and unique names
write.table(all_scores, "NES_scores_all_TFs.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Create a heatmap version
scores_matrix <- as.matrix(all_scores[,!(names(all_scores) %in% c("TF", "Original_TF"))])
rownames(scores_matrix) <- all_scores$TF

# Create breaks for color scale
max_abs <- max(abs(scores_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 100)

# Create annotation for duplicated TFs
row_annotation <- data.frame(
    Duplicated = factor(ifelse(all_scores$Original_TF %in% 
                              all_scores$Original_TF[duplicated(all_scores$Original_TF)],
                              "Yes", "No"))
)
rownames(row_annotation) <- all_scores$TF

ann_colors <- list(
    Duplicated = c(Yes = "red", No = "white")
)

# Create heatmap
library(pheatmap)
pheatmap(scores_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         filename = "NES_scores_heatmap.png",
         width = 12,
         height = 20)


