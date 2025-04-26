# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(enrichplot)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork)

# Running ORA on DGE results: HC vs. malignant 
# Using DGE filtered for OCR-Gene pairs 

# set working directory
setwd("ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR/")
load("coordination_pairs.RData")
# Sanity check
dim(MGUS_HC)  # 140
dim(MM_HC)    # 768

# Load DGE results

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
rna_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

OCR_gene_final <- rbind(MGUS_HC, MM_HC)
OCR_gene_final <- OCR_gene_final[!duplicated(OCR_gene_final$pairs),]
dim(OCR_gene_final) # 788

OCR_gene_final$rna_CyclindD1_MGUS <- 0
OCR_gene_final$rna_Mmset_MGUS <- 0
OCR_gene_final$rna_Trp53_MGUS <- 0
OCR_gene_final$rna_MIc_MGUS <- 0
OCR_gene_final$rna_CyclindD1_MM <- 0
OCR_gene_final$rna_Mmset_MM <- 0
OCR_gene_final$rna_Trp53_MM <- 0
OCR_gene_final$rna_MIc_MM <- 0
OCR_gene_final$rna_CyclindD1_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_CyclinD1_MGUS"]
OCR_gene_final$rna_Mmset_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Mmset_MGUS"]
OCR_gene_final$rna_Trp53_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Trp53_MGUS"]
OCR_gene_final$rna_MIc_MGUS[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_MIc_MGUS"]
OCR_gene_final$rna_CyclindD1_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_CyclinD1_MM"]
OCR_gene_final$rna_Mmset_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Mmset_MM"]
OCR_gene_final$rna_Trp53_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_Trp53_MM"]
OCR_gene_final$rna_MIc_MM[OCR_gene_final$gene_ID %in% rownames(dea_results)] <- dea_results[OCR_gene_final$gene_ID[OCR_gene_final$gene_ID %in% rownames(dea_results)], "sig_MIc_MM"]


# Split gene sets by direction
CyclinD1_MGUS_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_CyclindD1_MGUS == 1]
CyclinD1_MGUS_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_CyclindD1_MGUS == -1]
Mmset_MGUS_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Mmset_MGUS == 1]
Mmset_MGUS_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Mmset_MGUS == -1]
Trp53_MGUS_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Trp53_MGUS == 1]
Trp53_MGUS_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Trp53_MGUS == -1]
MIC_MGUS_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_MIc_MGUS == 1]
MIC_MGUS_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_MIc_MGUS == -1]
CyclinD1_MM_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_CyclindD1_MM == 1]
CyclinD1_MM_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_CyclindD1_MM == -1]
Mmset_MM_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Mmset_MM == 1]
Mmset_MM_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Mmset_MM == -1]
Trp53_MM_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Trp53_MM == 1]
Trp53_MM_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_Trp53_MM == -1]
MIC_MM_up <- OCR_gene_final$gene_ID[OCR_gene_final$rna_MIc_MM == 1]
MIC_MM_down <- OCR_gene_final$gene_ID[OCR_gene_final$rna_MIc_MM == -1]

# Load gene universe
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR")
universe <- readRDS(file = "DA_MGUS_OCR_gene.RDS")
universe <- unique(universe$gene_ID) # 1492


# Convert ENSEMBL to ENTREZ helper
convert_to_entrez <- function(ensembl_ids) {
  bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
    pull(ENTREZID) %>% na.omit() %>% unique()
}

universe_entrez <- convert_to_entrez(universe_ensembl)

# Create gene sets
gene_lists <- list(
CyclinD1_MGUS_up = CyclinD1_MGUS_up,
CyclinD1_MGUS_down = CyclinD1_MGUS_down,
Mmset_MGUS_up = Mmset_MGUS_up,
Mmset_MGUS_down = Mmset_MGUS_down,
Trp53_MGUS_up = Trp53_MGUS_up,  
Trp53_MGUS_down = Trp53_MGUS_down,
MIC_MGUS_up = MIC_MGUS_up,
MIC_MGUS_down = MIC_MGUS_down,
CyclinD1_MM_up = CyclinD1_MM_up,
CyclinD1_MM_down = CyclinD1_MM_down,
Mmset_MM_up = Mmset_MM_up,
Mmset_MM_down = Mmset_MM_down,
Trp53_MM_up = Trp53_MM_up,
Trp53_MM_down = Trp53_MM_down,
MIC_MM_up = MIC_MM_up,
MIC_MM_down = MIC_MM_down

)


# ORA function (with setReadable for KEGG + Hallmark)
run_ora <- function(ensembl_ids, name_prefix) {

  go_bp <- enrichGO(
    gene = ensembl_ids,
    universe = universe,
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE,
    minGSSize = 10,
    maxGSSize = 500
  )

  list(GO_BP = go_bp)
}

# Run ORA
ora_results <- list()
for (name in names(gene_lists)) {
  ensembl_ids <- gene_lists[[name]]
  ora_results[[name]] <- run_ora(ensembl_ids, name)
}


# Save Excel
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/Coordination/Gene_OCR/RNA/Enrichment/")
wb <- createWorkbook()
add_results_to_sheet <- function(wb, result_list, prefix) {
  for (category in names(result_list)) {
    df <- as.data.frame(result_list[[category]])
    if (nrow(df) > 0) {
      sheet_name <- substr(paste0(prefix, "_", category), 1, 31)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, df)
    }
  }
}
for (group in names(ora_results)) {
  add_results_to_sheet(wb, ora_results[[group]], group)
}
saveWorkbook(wb, file = "./ORA_Results_SPLIT_universeOCRGENES_OCRDA.xlsx", overwrite = TRUE)

# Plotting
pdf("ORA_Dot_Cnet_Plots_SPLIT_universeOCRGENES_OCRDA.pdf", width = 14, height = 7)
for (group in names(ora_results)) {
  for (category in names(ora_results[[group]])) {
    result <- ora_results[[group]][[category]]
    if (is.null(result) || nrow(result) == 0) next

    p1 <- tryCatch({
      dotplot(result, showCategory = 10, title = paste(group, category, "- Dotplot")) +
        theme(plot.title = element_text(hjust = 0.5))
    }, error = function(e) NULL)

    p2 <- tryCatch({
      cnetplot(result, showCategory = 5, circular = FALSE, colorEdge = TRUE) +
        ggtitle(paste(group, category, "- Cnetplot"))
    }, error = function(e) NULL)

    if (!is.null(p1) && !is.null(p2)) {
      print(p1 + p2)
    } else if (!is.null(p1)) {
      print(p1)
    } else if (!is.null(p2)) {
      print(p2)
    }
  }
}
dev.off()
