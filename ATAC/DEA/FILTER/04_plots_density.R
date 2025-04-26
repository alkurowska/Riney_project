#####################
### density plots ###
#####################

#####################
## Malignant vs HC ##
#####################


# Function
density_plot_2 <- function(toPlot, contrast) {
    p <- ggplot2::ggplot(toPlot, aes(x = MGUS_HC, fill = "MGUS_HC")) +
    geom_density(alpha = 0.9) +
    geom_density(aes(x = SMM_HC, fill = "SMM_HC"), alpha = 0.9) +
    geom_density(aes(x = MM_HC, fill = "MM_HC"), alpha = 0.9) +
    xlab("log2FC") + 
    ylab("Density") +
    # dash line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("MGUS_HC" = "#FFCC33", "SMM_HC" = "#CC3333", "MM_HC" = "#CC6633"), labels = c("MGUS vs HC","SMM vs HC", "MM vs HC")) +
    guides(fill = guide_legend(title = "Transition", override.aes = aes(label = ""))) +
    labs(color = "") +
    ggtitle("Density of log2FC") +
    theme_minimal()
    
    ggsave(p, filename = paste0(contrast, "_density_logFC.png"), width = 10, height = 10, dpi = 300)
    }


# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
dea_results <- read.table("atac_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/No_score")
dea_res_noscore <- read.table("atac_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER")
load("final_peaks.RData")


setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/FILTER/PLOTS")
# Global density plot of logFC
toPlot <- data.frame(dea_results$logFC_MGUS_HC, dea_results$logFC_SMM_HC, dea_res_noscore$logFC_MM_HC)
colnames(toPlot) <- c("MGUS_HC", "SMM_HC", "MM_HC")
rownames(toPlot) <- rownames(dea_results)

density_plot_2(toPlot, "global")

# Density plot of significant in MGUS vs HC
toPlot <- data.frame(dea_results[MGUS_HC_filter,]$logFC_MGUS_HC,
                        dea_results[MGUS_HC_filter,]$logFC_SMM_HC,
                        dea_res_noscore[MGUS_HC_filter,]$logFC_MM_HC)

colnames(toPlot) <- c("MGUS_HC", "SMM_HC", "MM_HC")
rownames(toPlot) <- MGUS_HC_filter

density_plot_2(toPlot, "MGUS_HC")

# Density plot of significant in SMM vs HC
toPlot <- data.frame(dea_results[SMM_HC_filter,]$logFC_MGUS_HC,
                        dea_results[SMM_HC_filter,]$logFC_SMM_HC,
                        dea_res_noscore[SMM_HC_filter,]$logFC_MM_HC)

colnames(toPlot) <- c("MGUS_HC", "SMM_HC", "MM_HC")
rownames(toPlot) <- SMM_HC_filter

density_plot_2(toPlot, "SMM_HC")

# Density plot of significant in MM vs HC
toPlot <- data.frame(dea_results[MM_HC_filter,]$logFC_MGUS_HC,
                        dea_results[MM_HC_filter,]$logFC_SMM_HC,
                        dea_res_noscore[MM_HC_filter,]$logFC_MM_HC)

colnames(toPlot) <- c("MGUS_HC", "SMM_HC", "MM_HC")
rownames(toPlot) <- MM_HC_filter

density_plot_2(toPlot, "MM_HC")