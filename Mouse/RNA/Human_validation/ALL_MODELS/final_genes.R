# Transalte mouse orthologs

# LOAD FINAL human genes
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")
final_genes <- unique(c(MGUS_HC, SMM_HC, MM_HC))

# LOAD all of the human genes
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/")
load("norm_to_voom.RData")
gene_anno_human <- gene_anno
gene_anno_human <- gene_anno_human[gene_anno_human$gene_id%in%rownames(d1_norm$counts),]
dim(gene_anno_human) # 24256 genes

all_genes <- gene_anno_human$gene_id

# Get mouse genes
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix/")
load("norm_to_voom.RData")
gene_anno_mouse <- read.table("gene_annotation.txt", sep = "\t", header = TRUE)
gene_anno_mouse <- gene_anno_mouse[gene_anno_mouse$ensembl_gene_id%in%rownames(d1_norm$counts),]
dim(gene_anno_mouse) # 19423 genes


# ortologs 

##### Convert human IDs to homolog mouse IDs with BIOMART #####
library(biomaRt)
library('org.Hs.eg.db') # Human
library('org.Mm.eg.db') # Mouse

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genes_all <- getLDS(attributes = "ensembl_gene_id",
       filters = "ensembl_gene_id", values = all_genes,
       mart = human, uniqueRows = T,
       attributesL = c("ensembl_gene_id"),
       martL = mouse)

colnames(genes_all) <- c("Human_ensembl_gene_id", "Mouse_ensembl_gene_id")
dim(genes_all) # 15586 genes

# Keep only the expressed mouse genes
genes_all <- genes_all[genes_all$Mouse_ensembl_gene_id %in% gene_anno_mouse$ensembl_gene_id,]
dim(genes_all) # 12701 genes

###>------------ it gives many duplicates

# Save homologs 
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/ALL_MODELS")
write.table(genes_all, file = "human_mouse_orthologs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Check final genes 
length(x = final_genes) # 1298 total genes 
final <- genes_all[genes_all$Human_ensembl_gene_id %in% final_genes,]
dim(final) # 978

length(unique(final$Human_ensembl_gene_id)) # 927 genes
length(unique(final$Mouse_ensembl_gene_id)) # 974 genes

# How much orthologs we have in general?
length(unique(genes_all$Human_ensembl_gene_id)) # 12295 genes
length(unique(genes_all$Mouse_ensembl_gene_id)) # 12353 genes

# How many genes are in the final genes list?
#927/1298 -> 0.7141757 71% of final genes
#12295/24256 -> 0.5068849 51% of all human genes

# Cor plots for each model across the 2 contrasts HC > MGUS and HC to MM
# Load the dea data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
human_noscore <- read.table("rna_dea_results_noscore.txt", sep = "\t", header = TRUE)
human_noscore <- human_noscore[unique(final$Human_ensembl_gene_id),]

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
human_gsva <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
human_gsva <- human_gsva[unique(final$Human_ensembl_gene_id),]

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/ALL_MODELS")
mouse_dea <- read.table("rna_dea_results.txt", sep = "\t", header = TRUE)
mouse_dea <- mouse_dea[unique(final$Mouse_ensembl_gene_id),]

# MGUS_HC conditions
final$mouse_MGUS_HC <- lapply(final$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "logFC_MGUS_HC"]
})

# MM_HC conditions
final$mouse_MM_HC <- lapply(final$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "logFC_MM_HC"]
})

# HUMAN
final$MGUS_HC <- lapply(final$Human_ensembl_gene_id, function(x) {
    human_gsva[x, "logFC_MGUS_HC"]
})
final$MM_HC <- lapply(final$Human_ensembl_gene_id, function(x) {
    human_noscore[x, "logFC_MM_HC"]
})


# Convert all list columns to vectors
final <- data.frame(lapply(final, unlist))

# For two contrasts: MGUS_HC and MM_HC make correlation of logFC between mouse and human homologs for each model separately
# MGUS_HC
MGUS_HC <- final[,c("MGUS_HC", "mouse_MGUS_HC")]
MM_HC <- final[,c("MM_HC", "mouse_MM_HC")]
log_list <- list(
    MGUS_HC = MGUS_HC,
    MM_HC = MM_HC
)

# Function
cor_plot <- function(toPlot, contrast) {
    data <- toPlot
    rownames(data) <- 1:nrow(data)
    colnames(data) <- c("Human", "Mouse")

    correlation_coefficient <- cor(data$Human, data$Mouse, method = "spearman")

    p <-  ggplot(data, aes(x = Human, y = Mouse)) +
            geom_point() +  # Plot points
            # add line in x = 0 
            geom_vline(xintercept = 0, color = "black") +
            geom_hline(yintercept = 0, color = "black") +
            annotate("text", x = -4.5, y = 4.5, label = paste("rho=", round(correlation_coefficient, 5))) +
            xlab("Mouse") + 
            ylab("Human") +
            xlim(-5, 5) +
            ylim(-5, 5) +
            ggtitle("Correlation of log2FC") +
            theme_minimal()
    ggsave(p, filename = paste0("correlation_",contrast,".png"), width = 10, height = 10, dpi = 300)
    }

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/ALL_MODELS/")
cor_plot(MGUS_HC, "MGUS_HC")
cor_plot(MM_HC, "MM_HC")


# For each model plot sankey plot separately for up and down genes

# GET ALL UP AND DOWN HUAMN FINAL GENES with mouse homologs
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

MGUS_HC_up <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==1,])]
SMM_HC_up <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==1,])]
MM_HC_up <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==1,])]

MGUS_HC_down <- MGUS_HC[MGUS_HC%in%rownames(dea_results[dea_results$sig_MGUS_HC==-1,])]
SMM_HC_down <- SMM_HC[SMM_HC%in%rownames(dea_results[dea_results$sig_SMM_HC==-1,])]
MM_HC_down <- MM_HC[MM_HC%in%rownames(dea_res_noscore[dea_res_noscore$sig_MM_HC==-1,])]

# UP in all
all_up <- unique(c(MGUS_HC_up, SMM_HC_up, MM_HC_up))

# Down in all
all_down <- unique(c(MGUS_HC_down, SMM_HC_down, MM_HC_down))

# Keep only the final genes with homologs 
final_up <- final[final$Human_ensembl_gene_id %in% all_up,]
final_down <- final[final$Human_ensembl_gene_id %in% all_down,]

length(unique(final_up$Human_ensembl_gene_id)) # 254 genes
length(unique(final_down$Human_ensembl_gene_id)) # 663 genes

# Keep only the unique mouse genes
final_up <- final_up[!duplicated(final_up$Mouse_ensembl_gene_id),]
final_down <- final_down[!duplicated(final_down$Mouse_ensembl_gene_id),]


# MGUS_HC conditions
final_up$CyclinD1_MGUS_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_CyclinD1_MGUS"]
})

final_up$MIc_MGUS_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_MIc_MGUS"]
})

final_up$Trp53_MGUS_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Trp53_MGUS"]
})

final_up$Mmset_MGUS_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Mmset_MGUS"]
})

# MM_HC conditions
final_up$CyclinD1_MM_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_CyclinD1_MM"]
})

final_up$MIc_MM_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_MIc_MM"]
})

final_up$Trp53_MM_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Trp53_MM"]
})

final_up$Mmset_MM_HC <- lapply(final_up$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Mmset_MM"]
})

final_up$MGUS_HC <- lapply(final_up$Human_ensembl_gene_id, function(x) {
    human_gsva[x, "sig_MGUS_HC"]
})
final_up$MM_HC <- lapply(final_up$Human_ensembl_gene_id, function(x) {
    human_noscore[x, "sig_MM_HC"]
})

# Convert all list columns to vectors
final_up <- data.frame(lapply(final_up, unlist))


# Sanky plot
library(ggplot2)
library(ggalluvial)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
#### PREP THE DATA ####

models <- c("CyclinD1", "MIc", "Trp53", "Mmset")
for(i in 1:length(models)) {
    model <- models[i]
    toPlot <- final_up[,c(paste0(model, "_MGUS_HC"), paste0(model, "_MM_HC"))]
    colnames(toPlot) <- c("MGUS_vs_HC", "MM_vs_HC")
    # Change 0s to "non differential"
    toPlot[toPlot == 0] <- "Non-differential"
    #Change -1 to "Down"
    toPlot[toPlot == -1] <- "Down-regulated"
    #Change 1 to "Up"
    toPlot[toPlot == 1] <- "Up-regulated"

    toPlot$Gene <- final_up$Mouse_ensembl_gene_id
    # Reshape for ggalluvial
    genes_long <- tidyr::pivot_longer(
    toPlot,
    cols = c(MGUS_vs_HC, MM_vs_HC),
    names_to = "Condition",
    values_to = "Category")

    genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "MM_vs_HC"))

    # Plot
    p <- ggplot(genes_long, aes(
        x = Condition, stratum = Category, alluvium = Gene,
        fill = Category, label = Category
        )) +
        geom_flow(stat = "alluvium") +
        geom_stratum() +
        theme_classic() +
    labs(x = NULL, y = "Gene Count") +
    scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
    ggsave(p, filename = paste0(model, "_sankey_up.png"), width = 5, height = 5, dpi = 300)

}

# DOWN

# MGUS_HC conditions
final_down$CyclinD1_MGUS_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_CyclinD1_MGUS"]
})

final_down$MIc_MGUS_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_MIc_MGUS"]
})

final_down$Trp53_MGUS_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Trp53_MGUS"]
})

final_down$Mmset_MGUS_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Mmset_MGUS"]
})

# MM_HC conditions
final_down$CyclinD1_MM_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_CyclinD1_MM"]
})

final_down$MIc_MM_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_MIc_MM"]
})

final_down$Trp53_MM_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Trp53_MM"]
})

final_down$Mmset_MM_HC <- lapply(final_down$Mouse_ensembl_gene_id, function(x) {
    mouse_dea[x, "sig_Mmset_MM"]
})

final_down$MGUS_HC <- lapply(final_down$Human_ensembl_gene_id, function(x) {
    human_gsva[x, "sig_MGUS_HC"]
})
final_down$MM_HC <- lapply(final_down$Human_ensembl_gene_id, function(x) {
    human_noscore[x, "sig_MM_HC"]
})

# Convert all list columns to vectors
final_down <- data.frame(lapply(final_down, unlist))


# Sanky plot
library(ggplot2)
library(ggalluvial)
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
#### PREP THE DATA ####

models <- c("CyclinD1", "MIc", "Trp53", "Mmset")
for(i in 1:length(models)) {
    model <- models[i]
    toPlot <- final_up[,c(paste0(model, "_MGUS_HC"), paste0(model, "_MM_HC"))]
    colnames(toPlot) <- c("MGUS_vs_HC", "MM_vs_HC")
    # Change 0s to "non differential"
    toPlot[toPlot == 0] <- "Non-differential"
    #Change -1 to "Down"
    toPlot[toPlot == -1] <- "Down-regulated"
    #Change 1 to "Up"
    toPlot[toPlot == 1] <- "Up-regulated"

    toPlot$Gene <- final_up$Mouse_ensembl_gene_id
    # Reshape for ggalluvial
    genes_long <- tidyr::pivot_longer(
    toPlot,
    cols = c(MGUS_vs_HC, MM_vs_HC),
    names_to = "Condition",
    values_to = "Category")

    genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "MM_vs_HC"))

    # Plot
    p <- ggplot(genes_long, aes(
        x = Condition, stratum = Category, alluvium = Gene,
        fill = Category, label = Category
        )) +
        geom_flow(stat = "alluvium") +
        geom_stratum() +
        theme_classic() +
    labs(x = NULL, y = "Gene Count") +
    scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
    ggsave(p, filename = paste0(model, "_sankey_down.png"), width = 5, height = 5, dpi = 300)

}

