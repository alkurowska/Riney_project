###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

library(UpSetR)

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

CyclinD1_MGUS_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS_HC <- rownames(dea_results[dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS_HC <- rownames(dea_results[dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS_HC <- rownames(dea_results[dea_results$sig_MIc_MGUS != 0,])
CyclinD1_MM_HC <- rownames(dea_results[dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM_HC <- rownames(dea_results[dea_results$sig_Mmset_MM != 0,])
Trp53_MM_HC <- rownames(dea_results[dea_results$sig_Trp53_MM != 0,])
MIc_MM_HC <- rownames(dea_results[dea_results$sig_MIc_MM != 0,])

# Directions 

CyclinD1_MGUS_HC_up <- CyclinD1_MGUS_HC[CyclinD1_MGUS_HC%in%rownames(dea_results[dea_results$sig_CyclinD1_MGUS==1,])]
CyclinD1_MM_HC_up <- CyclinD1_MM_HC[CyclinD1_MM_HC%in%rownames(dea_results[dea_results$sig_CyclinD1_MM==1,])]
Mmset_MGUS_HC_up <- Mmset_MGUS_HC[Mmset_MGUS_HC%in%rownames(dea_results[dea_results$sig_Mmset_MGUS==1,])]
Mmset_MM_HC_up <- Mmset_MM_HC[Mmset_MM_HC%in%rownames(dea_results[dea_results$sig_Mmset_MM==1,])]
Trp53_MGUS_HC_up <- Trp53_MGUS_HC[Trp53_MGUS_HC%in%rownames(dea_results[dea_results$sig_Trp53_MGUS==1,])]
Trp53_MM_HC_up <- Trp53_MM_HC[Trp53_MM_HC%in%rownames(dea_results[dea_results$sig_Trp53_MM==1,])]
MIc_MGUS_HC_up <- MIc_MGUS_HC[MIc_MGUS_HC%in%rownames(dea_results[dea_results$sig_MIc_MGUS==1,])]
MIc_MM_HC_up <- MIc_MM_HC[MIc_MM_HC%in%rownames(dea_results[dea_results$sig_MIc_MM==1,])]

CyclinD1_MGUS_HC_down <- CyclinD1_MGUS_HC[CyclinD1_MGUS_HC%in%rownames(dea_results[dea_results$sig_CyclinD1_MGUS==-1,])]
CyclinD1_MM_HC_down <- CyclinD1_MM_HC[CyclinD1_MM_HC%in%rownames(dea_results[dea_results$sig_CyclinD1_MM==-1,])]
Mmset_MGUS_HC_down <- Mmset_MGUS_HC[Mmset_MGUS_HC%in%rownames(dea_results[dea_results$sig_Mmset_MGUS==-1,])]
Mmset_MM_HC_down <- Mmset_MM_HC[Mmset_MM_HC%in%rownames(dea_results[dea_results$sig_Mmset_MM==-1,])]
Trp53_MGUS_HC_down <- Trp53_MGUS_HC[Trp53_MGUS_HC%in%rownames(dea_results[dea_results$sig_Trp53_MGUS==-1,])]
Trp53_MM_HC_down <- Trp53_MM_HC[Trp53_MM_HC%in%rownames(dea_results[dea_results$sig_Trp53_MM==-1,])]
MIc_MGUS_HC_down <- MIc_MGUS_HC[MIc_MGUS_HC%in%rownames(dea_results[dea_results$sig_MIc_MGUS==-1,])]
MIc_MM_HC_down <- MIc_MM_HC[MIc_MM_HC%in%rownames(dea_results[dea_results$sig_MIc_MM==-1,])]

sets <- unique(c(CyclinD1_MGUS_HC, CyclinD1_MM_HC, Mmset_MGUS_HC, Mmset_MM_HC, Trp53_MGUS_HC, Trp53_MM_HC, MIc_MGUS_HC, MIc_MM_HC))

# Create a binary matrix of membership
all_genes <- unique(unlist(sets))

membership_matrix <- data.frame(
  Gene = all_genes,
  CyclinD1_MGUS_HC_up = all_genes %in% CyclinD1_MGUS_HC_up,
  CyclinD1_MM_HC_up = all_genes %in% CyclinD1_MM_HC_up,
  Mmset_MGUS_HC_up = all_genes %in% Mmset_MGUS_HC_up,
  Mmset_MM_HC_up = all_genes %in% Mmset_MM_HC_up,
  Trp53_MGUS_HC_up = all_genes %in% Trp53_MGUS_HC_up,
  Trp53_MM_HC_up = all_genes %in% Trp53_MM_HC_up,
  MIc_MGUS_HC_up = all_genes %in% MIc_MGUS_HC_up,
  MIc_MM_HC_up = all_genes %in% MIc_MM_HC_up
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]
colnames(upset_data) <- c("CyclinD1 MGUS vs. HC up", "CyclinD1 MM vs. HC up", "Mmset MGUS vs HC up",
                            "Mmset MM vs HC up", "Trp53 MGUS vs HC up", "Trp53 MM vs HC up",
                            "MIc MGUS vs HC up", "MIc MM vs HC up")
# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_directions_up.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("CyclinD1 MGUS vs. HC up", "CyclinD1 MM vs. HC up", "Mmset MGUS vs HC up",
                            "Mmset MM vs HC up", "Trp53 MGUS vs HC up", "Trp53 MM vs HC up",
                            "MIc MGUS vs HC up", "MIc MM vs HC up"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

# down 
membership_matrix_down <- data.frame(
  Gene = all_genes,
  CyclinD1_MGUS_HC_down = all_genes %in% CyclinD1_MGUS_HC_down,
  CyclinD1_MM_HC_down = all_genes %in% CyclinD1_MM_HC_down,
  Mmset_MGUS_HC_down = all_genes %in% Mmset_MGUS_HC_down,
  Mmset_MM_HC_down = all_genes %in% Mmset_MM_HC_down,
  Trp53_MGUS_HC_down = all_genes %in% Trp53_MGUS_HC_down,
  Trp53_MM_HC_down = all_genes %in% Trp53_MM_HC_down,
  MIc_MGUS_HC_down = all_genes %in% MIc_MGUS_HC_down,
  MIc_MM_HC_down = all_genes %in% MIc_MM_HC_down
)

# Convert logical values to numeric for plotting
membership_matrix_down[, -1] <- lapply(membership_matrix_down[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data_down <- membership_matrix_down[, -1]

colnames(upset_data_down) <- c("CyclinD1 MGUS vs. HC down", "CyclinD1 MM vs. HC down", "Mmset MGUS vs HC down",
                            "Mmset MM vs HC down", "Trp53 MGUS vs HC down", "Trp53 MM vs HC down",
                            "MIc MGUS vs HC down", "MIc MM vs HC down")

# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_directions_down.png", width=2500, height=2000, res=300)
upset(upset_data_down, sets = c("CyclinD1 MGUS vs. HC down", "CyclinD1 MM vs. HC down", "Mmset MGUS vs HC down",
                            "Mmset MM vs HC down", "Trp53 MGUS vs HC down", "Trp53 MM vs HC down",
                            "MIc MGUS vs HC down", "MIc MM vs HC down"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

# General
membership_matrix_general <- data.frame(
  Gene = all_genes,
  CyclinD1_MGUS_HC = all_genes %in% CyclinD1_MGUS_HC,
  CyclinD1_MM_HC = all_genes %in% CyclinD1_MM_HC,
  Mmset_MGUS_HC = all_genes %in% Mmset_MGUS_HC,
  Mmset_MM_HC = all_genes %in% Mmset_MM_HC,
  Trp53_MGUS_HC = all_genes %in% Trp53_MGUS_HC,
  Trp53_MM_HC = all_genes %in% Trp53_MM_HC,
  MIc_MGUS_HC = all_genes %in% MIc_MGUS_HC,
  MIc_MM_HC = all_genes %in% MIc_MM_HC
)

# Convert logical values to numeric for plotting
membership_matrix_general[, -1] <- lapply(membership_matrix_general[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data_general <- membership_matrix_general[, -1]
colnames(upset_data_general) <- c("CyclinD1 MGUS vs. HC", "CyclinD1 MM vs. HC", "Mmset MGUS vs HC",
                            "Mmset MM vs HC", "Trp53 MGUS vs HC", "Trp53 MM vs HC",
                            "MIc MGUS vs HC", "MIc MM vs HC")

# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_directions_general.png", width=2500, height=2000, res=300)
upset(upset_data_general, sets = c("CyclinD1 MGUS vs. HC", "CyclinD1 MM vs. HC", "Mmset MGUS vs HC",
                            "Mmset MM vs HC", "Trp53 MGUS vs HC", "Trp53 MM vs HC",
                            "MIc MGUS vs HC", "MIc MM vs HC"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()


# plot without dividing per contrast
CyclinD1 <- unique(c(CyclinD1_MGUS_HC, CyclinD1_MM_HC))
Mmset <- unique(c(Mmset_MGUS_HC, Mmset_MM_HC))
Trp53 <- unique(c(Trp53_MGUS_HC, Trp53_MM_HC))
MIc <- unique(c(MIc_MGUS_HC, MIc_MM_HC))

sets <- unique(c(CyclinD1, Mmset, Trp53, MIc))
# Create a binary matrix of membership
all_genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = all_genes,
  CyclinD1 = all_genes %in% CyclinD1,
  Mmset = all_genes %in% Mmset,
  Trp53 = all_genes %in% Trp53,
  MIc = all_genes %in% MIc
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("CyclinD1", "Mmset", "Trp53", "MIc")


# Generate the UpSet plot

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_directions_general_no_contrast.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("CyclinD1", "Mmset", "Trp53", "MIc"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

# Separate for MGUS without directions and MM without directions
CyclinD1_MGUS <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS != 0,])
Mmset_MGUS <- rownames(dea_results[dea_results$sig_Mmset_MGUS != 0,])
Trp53_MGUS <- rownames(dea_results[dea_results$sig_Trp53_MGUS != 0,])
MIc_MGUS <- rownames(dea_results[dea_results$sig_MIc_MGUS != 0,])

CyclinD1_MM <- rownames(dea_results[dea_results$sig_CyclinD1_MM != 0,])
Mmset_MM <- rownames(dea_results[dea_results$sig_Mmset_MM != 0,])
Trp53_MM <- rownames(dea_results[dea_results$sig_Trp53_MM != 0,])
MIc_MM <- rownames(dea_results[dea_results$sig_MIc_MM != 0,])

sets <- unique(c(CyclinD1_MGUS, CyclinD1_MM, Mmset_MGUS, Mmset_MM, Trp53_MGUS, Trp53_MM, MIc_MGUS, MIc_MM))
# Create a binary matrix of membership
all_genes <- unique(unlist(sets))

membership_matrix <- data.frame(
  Gene = all_genes,
  CyclinD1_MGUS = all_genes %in% CyclinD1_MGUS,
  Mmset_MGUS = all_genes %in% Mmset_MGUS,
  Trp53_MGUS = all_genes %in% Trp53_MGUS,
  MIc_MGUS = all_genes %in% MIc_MGUS
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column


upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("CyclinD1 MGUS", "Mmset MGUS", "Trp53 MGUS", "MIc MGUS")
# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_MGUS.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("CyclinD1 MGUS", "Mmset MGUS", "Trp53 MGUS", "MIc MGUS"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

# Separate for MGUS without directions and MM without directions
sets <- unique(c(CyclinD1_MM, Mmset_MM, Trp53_MM, MIc_MM))
# Create a binary matrix of membership
all_genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = all_genes,
  CyclinD1_MM = all_genes %in% CyclinD1_MM,
  Mmset_MM = all_genes %in% Mmset_MM,
  Trp53_MM = all_genes %in% Trp53_MM,
  MIc_MM = all_genes %in% MIc_MM
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column

upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("CyclinD1 MM", "Mmset MM", "Trp53 MM", "MIc MM")
# Generate the UpSet plot
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
png("upset_plot_MM.png", width=2500, height=2000, res=300)
upset(upset_data, sets = c("CyclinD1 MM", "Mmset MM", "Trp53 MM", "MIc MM"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()

