###########################
### GLOBAL SET of GENES ###
###########################

## model 1 - no score
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

MM_HC_no <- rownames(dea_res_noscore)[which(dea_res_noscore$sig_MM_HC != 0)]
SMM_HC_no <- rownames(dea_res_noscore)[which(dea_res_noscore$sig_SMM_HC != 0)]
MGUS_HC_no <- rownames(dea_res_noscore)[which(dea_res_noscore$sig_MGUS_HC != 0)]


global <- unique(c(MM_HC_no, SMM_HC_no, MGUS_HC_no)) #5884
length(global)

###########################

## model 2 - with score
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

MM_HC <- rownames(dea_results)[which(dea_results$sig_MM_HC != 0)]
SMM_HC <- rownames(dea_results)[which(dea_results$sig_SMM_HC != 0)]
MGUS_HC <- rownames(dea_results)[which(dea_results$sig_MGUS_HC != 0)]


# FILTER each contrast as an intersection with the global set
MM_HC_no_filter <- MM_HC_no[MM_HC_no %in% global] #5068
SMM_HC_filter <- SMM_HC[SMM_HC %in% global] #2431
MGUS_HC_filter <- MGUS_HC[MGUS_HC %in% global] #1263


# Venn diagram
library(ggVennDiagram)
library(ggplot2)


# Venn diagram
# Create a Named List of Sets
venn_data <- list(MM_HC_no = MM_HC_no, GLOBAL = global)

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
# Plot the Venn Diagram with Custom Colors
png("venn_plot_MM.png", width=1500, height=1000, res=300)
ggVennDiagram(venn_data, label_alpha = 0) +  
  # Make the plot horizontal
  coord_flip() +
  scale_fill_gradient(low = "grey", high = "purple",  
                      limits = c(0, length(global))) +  # Custom colors
  theme()
dev.off()


# Create a Named List of Sets
venn_data <- list(SMM_HC = SMM_HC, GLOBAL = global)

# Plot the Venn Diagram with Custom Colors
png("venn_plot_SMM.png", width=1500, height=1000, res=300)
ggVennDiagram(venn_data, label_alpha = 0) +  
  # Make the plot horizontal
  coord_flip() +
  scale_fill_gradient(low = "grey", high = "purple",  
                      limits = c(0, length(global))) +  # Custom colors
  theme()
dev.off()



# Create a Named List of Sets
venn_data <- list(MGUS_HC = MGUS_HC, GLOBAL = global)

# Plot the Venn Diagram with Custom Colors
png("venn_plot_MGUS.png", width=1500, height=1000, res=300)
ggVennDiagram(venn_data, label_alpha = 0) +  
  # Make the plot horizontal
  coord_flip() +
  scale_fill_gradient(low = "grey", high = "purple",  
                      limits = c(0, length(global))) +  # Custom colors
  theme()
dev.off()


## Venn diagram of 3 final sets
venn_data <- list(MM_HC = MM_HC_no_filter, SMM_HC = SMM_HC_filter, MGUS_HC = MGUS_HC_filter)

# Plot the Venn Diagram with Custom Colors
png("venn_plot_final.png", width=2000, height=2000, res=300)
ggVennDiagram(venn_data, label_alpha = 0) +  
  # Make the plot horizontal
  scale_fill_gradient(low = "grey", high = "purple",  
                      limits = c(0, length(MM_HC_no_filter))) +  # Custom colors
  theme()
dev.off()


# Save the data
MM_HC_filter <-  MM_HC_no_filter
save(MM_HC_filter, SMM_HC_filter, MGUS_HC_filter, global, file="final_genes.RData")