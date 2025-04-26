# Read homologs 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
orthologs <- read.table("human_mouse_orthologs.txt", sep = "\t", header = TRUE)

setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")

# Cor plots for each model across the 2 contrasts HC > MGUS and HC to MM
# Load the dea data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
human_noscore <- read.table("rna_dea_results_noscore.txt", sep = "\t", header = TRUE)
human_noscore <- human_noscore[unique(orthologs$Human_ensembl_gene_id),]

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
human_gsva <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
human_gsva <- human_gsva[unique(orthologs$Human_ensembl_gene_id),]

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA")
mouse_dea <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
mouse_dea <- mouse_dea[unique(orthologs$Mouse_ensembl_gene_id),]

orthologs$CyclinD1_MGUS_HC <- NA
orthologs$MIc_MGUS_HC <- NA
orthologs$Trp53_MGUS_HC <- NA
orthologs$Mmset_MGUS_HC <- NA
orthologs$CyclinD1_MM_HC <- NA
orthologs$MIc_MM_HC <- NA
orthologs$Trp53_MM_HC <- NA
orthologs$Mmset_MM_HC <- NA

for(i in 1:nrow(orthologs)) {
    orthologs$CyclinD1_MGUS_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_CyclinD1_MGUS"]
    orthologs$MIc_MGUS_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_MIc_MGUS"]
    orthologs$Trp53_MGUS_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_Trp53_MGUS"]
    orthologs$Mmset_MGUS_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_Mmset_MGUS"]
    orthologs$CyclinD1_MM_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_CyclinD1_MM"]
    orthologs$MIc_MM_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_MIc_MM"]
    orthologs$Trp53_MM_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_Trp53_MM"]
    orthologs$Mmset_MM_HC[i] <- mouse_dea[rownames(mouse_dea)==orthologs$Mouse_ensembl_gene_id[i], "logFC_Mmset_MM"]
}

orthologs$MM_HC <- NA
orthologs$SMM_HC <- NA
orthologs$MGUS_HC <- NA

for(i in 1:nrow(orthologs)){
    orthologs$MGUS_HC[i] <- human_gsva[rownames(human_gsva)==orthologs$Human_ensembl_gene_id[i], "logFC_MGUS_HC"]
    orthologs$MM_HC[i] <- human_noscore[rownames(human_noscore)==orthologs$Human_ensembl_gene_id[i], "logFC_MM_HC"]
    orthologs$SMM_HC[i] <- human_gsva[rownames(human_gsva)==orthologs$Human_ensembl_gene_id[i], "logFC_SMM_HC"]
}

# For two contrasts: MGUS_HC and MM_HC make correlation of logFC between mouse and human homologs for each model separately
# MGUS_HC # 137
MGUSvsHC <- orthologs[orthologs$Human_ensembl_gene_id%in%MGUS_HC,c("MGUS_HC", "CyclinD1_MGUS_HC", "MIc_MGUS_HC", "Trp53_MGUS_HC", "Mmset_MGUS_HC")]
dim(SMMvsHC) #338
SMMvsHC <- orthologs[orthologs$Human_ensembl_gene_id%in%SMM_HC,c("SMM_HC", "CyclinD1_MM_HC", "MIc_MM_HC", "Trp53_MM_HC", "Mmset_MM_HC")]
# MM_HC # 944
MMvsHC <- orthologs[orthologs$Human_ensembl_gene_id%in%MM_HC,c("MM_HC", "CyclinD1_MM_HC", "MIc_MM_HC", "Trp53_MM_HC", "Mmset_MM_HC")]
log_list <- list(
    MGUS_HC = MGUSvsHC,
    MM_HC = MMvsHC,
    SMM_HC = SMMvsHC
)
#SMM_HC 
table(unique(orthologs$Human_ensembl_gene_id)%in%MGUS_HC)

library(ggplot2)
# Function
cor_plot <- function(toPlot, contrast, model) {
    data <- toPlot
    rownames(data) <- 1:nrow(data)
    colnames(data) <- c("Human", "Mouse")

    correlation_coefficient <- cor(data$Human, data$Mouse, method = "spearman")

    p <-  ggplot(data, aes(x = Human, y = Mouse)) +
            geom_point() +  # Plot points
            # add line in x = 0 
            geom_vline(xintercept = 0, color = "black") +
            geom_hline(yintercept = 0, color = "black") +
            annotate("text", x = -4, y = 4.5, label = paste("rho=", round(correlation_coefficient, 5))) +
            xlab("Human") + 
            ylab("Mouse") +
            xlim(-5, 5) +
            ylim(-5, 5) +
            ggtitle(paste0("Correlation of log2FC - ", model)) +
            theme_minimal()
    ggsave(p, filename = paste0("correlation_",model,"_",contrast,".png"), width = 5, height = 5, dpi = 300)
    }

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation")
for(i in 1:length(log_list)) {
    contrast <- names(log_list)[i]
    for(j in 1:(length(log_list[[i]])-1)) {
        model <- gsub("_.*","",names(log_list[[i]])[1+j])
        toPlot <- cbind(log_list[[i]][1], log_list[[i]][1+j])
        cor_plot(toPlot, contrast, model)
    }
}
