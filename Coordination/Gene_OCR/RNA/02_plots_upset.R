###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

library(UpSetR)

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")

# Generate the UpSet plot

sets <- list(MGUS_HC, SMM_HC, MM_HC)

# Create a binary matrix of membership
all_genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = all_genes,
  MGUS_HC = all_genes %in% MGUS_HC,
  SMM_HC = all_genes %in% SMM_HC,
  MM_HC = all_genes %in% MM_HC
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]
# Generate the UpSet plot

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")
png("upset_plot_general.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("MGUS_HC", "SMM_HC", "MM_HC"),
      order.by = "freq", 
      main.bar.color = "steelblue")

dev.off()


# Directions 

MGUS_HC_up <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==1,])]
SMM_HC_up <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==1,])]
MM_HC_up <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])]

MGUS_HC_down <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==-1,])]
SMM_HC_down <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==-1,])]
MM_HC_down <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])]

# Create a binary matrix of membership
all_genes <- unique(unlist(sets))

membership_matrix <- data.frame(
  Gene = all_genes,
  MGUS_HC_up = all_genes %in% MGUS_HC_up,
  SMM_HC_up = all_genes %in% SMM_HC_up,
  MM_HC_up = all_genes %in% MM_HC_up,
  MGUS_HC_down = all_genes %in% MGUS_HC_down,
  SMM_HC_down = all_genes %in% SMM_HC_down,
  MM_HC_down = all_genes %in% MM_HC_down
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]

# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")
png("upset_plot_directions.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("MGUS_HC_up", "SMM_HC_up", "MM_HC_up", "MGUS_HC_down", "SMM_HC_down", "MM_HC_down"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()