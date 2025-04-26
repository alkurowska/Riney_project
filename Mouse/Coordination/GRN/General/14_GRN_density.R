# Read GRN 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/GRN/General/promoters")
grn_prom <- read.table("OCR_gene_association.txt", header = T, sep = "\t")
# remove NA rows
grn_prom <- grn_prom[!is.na(grn_prom$peak),] # 11,138

toPlot <- grn_prom$spearman_rho

library(ggplot2)
library(reshape2)
toPlot_melted <- melt(toPlot)

# Create the density plot

p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot.png", p, width = 6, height = 4, dpi = 300)

# Density plot of significant 
# Get the significant 
sig <- grn_prom[grn_prom$spearman_adj.pval < 0.01,]

toPlot <- sig$spearman_rho

toPlot_melted <- melt(toPlot)

# Create the density plot
p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho of significant genes",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot_sig.png", p, width = 6, height = 4, dpi = 300)



# Repeat for enhancers
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/GRN/General/enhancers")
grn_enh <- read.table("OCR_gene_association.txt", header = T, sep = "\t")
# remove NA rows
grn_enh <- grn_enh[!is.na(grn_enh$peak),] 

toPlot <- grn_enh$spearman_rho

toPlot_melted <- melt(toPlot)

# Create the density plot
p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot_enh.png", p, width = 6, height = 4, dpi = 300)

# Density plot of significant
# Get the significant

sig <- grn_enh[grn_enh$spearman_adj.pval < 0.01,]

toPlot <- sig$spearman_rho

toPlot_melted <- melt(toPlot)

# Create the density plot

p <- ggplot(toPlot_melted, aes(x = value)) +
    geom_density(alpha = 0.5, color = "black", fill = "green") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
    labs(title = "Density plot of rho of significant genes",
       x = "Spearman's rho",
       y = "Density",
       fill = NULL)

# Save the plot
ggsave("rho_density_plot_sig_enh.png", p, width = 6, height = 4, dpi = 300)