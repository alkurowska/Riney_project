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
# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

CyclinD1_MGUS_HC_up <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS == 1,])
CyclinD1_MGUS_HC_down <- rownames(dea_results[dea_results$sig_CyclinD1_MGUS == -1,])
Mmset_MGUS_HC_up <- rownames(dea_results[dea_results$sig_Mmset_MGUS == 1,])
Mmset_MGUS_HC_down <- rownames(dea_results[dea_results$sig_Mmset_MGUS == -1,])
Trp53_MGUS_HC_up <- rownames(dea_results[dea_results$sig_Trp53_MGUS == 1,])
Trp53_MGUS_HC_down <- rownames(dea_results[dea_results$sig_Trp53_MGUS == -1,])
MIc_MGUS_HC_up <- rownames(dea_results[dea_results$sig_MIc_MGUS == 1,])
MIc_MGUS_HC_down <- rownames(dea_results[dea_results$sig_MIc_MGUS == -1,])
CyclinD1_MM_HC_up <- rownames(dea_results[dea_results$sig_CyclinD1_MM == 1,])
CyclinD1_MM_HC_down <- rownames(dea_results[dea_results$sig_CyclinD1_MM == -1,])
Mmset_MM_HC_up <- rownames(dea_results[dea_results$sig_Mmset_MM == 1,])
Mmset_MM_HC_down <- rownames(dea_results[dea_results$sig_Mmset_MM == -1,])
Trp53_MM_HC_up <- rownames(dea_results[dea_results$sig_Trp53_MM == 1,])
Trp53_MM_HC_down <- rownames(dea_results[dea_results$sig_Trp53_MM == -1,])
MIc_MM_HC_up <- rownames(dea_results[dea_results$sig_MIc_MM == 1,])
MIc_MM_HC_down <- rownames(dea_results[dea_results$sig_MIc_MM == -1,])


# Load gene universe - all differential
universe <- unique(c(
  CyclinD1_MGUS_HC_up, CyclinD1_MGUS_HC_down,
  Mmset_MGUS_HC_up, Mmset_MGUS_HC_down,
  Trp53_MGUS_HC_up, Trp53_MGUS_HC_down,
  MIc_MGUS_HC_up, MIc_MGUS_HC_down,
  CyclinD1_MM_HC_up, CyclinD1_MM_HC_down,
  Mmset_MM_HC_up, Mmset_MM_HC_down,
  Trp53_MM_HC_up, Trp53_MM_HC_down,
  MIc_MM_HC_up, MIc_MM_HC_down
))



# Load the HUMAN DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
human_dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
human_dea_results_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load overlaping orthologs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
orthologs_85 <- read.table("orthologs_85.txt", sep = "\t", header = TRUE)


# ADD dynamics 
# Add differential results
orthologs_85 <- merge(orthologs_85, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
orthologs_85 <- merge(orthologs_85, human_dea_results[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
orthologs_85$sig_MM_HC <- human_dea_results_noscore[orthologs_85$Human_ensembl_gene_id,]$sig_MM_HC

# Keep rows only with 0 and 1
orthologs_85_up <- orthologs_85[orthologs_85$sig_CyclinD1_MGUS %in% c(0, 1) & orthologs_85$sig_CyclinD1_MM %in% c(0, 1) & orthologs_85$sig_MIc_MM %in% c(0, 1) & orthologs_85$sig_MIc_MGUS %in% c(0, 1) & orthologs_85$sig_Mmset_MM %in% c(0, 1) & orthologs_85$sig_Mmset_MGUS %in% c(0, 1) & orthologs_85$sig_Trp53_MM %in% c(0, 1) & orthologs_85$sig_Trp53_MGUS %in% c(0, 1) & orthologs_85$sig_MGUS_HC %in% c(0, 1) & orthologs_85$sig_SMM_HC %in% c(0, 1) & orthologs_85$sig_MM_HC %in% c(0, 1),]
# remove rows with 0 for all human contrasts
orthologs_85_up <- orthologs_85_up[!(orthologs_85_up$sig_MGUS_HC == 0 & orthologs_85_up$sig_SMM_HC == 0 & orthologs_85_up$sig_MM_HC == 0),]

# Keep rows only with 0 and -1
orthologs_85_down <- orthologs_85[orthologs_85$sig_CyclinD1_MGUS %in% c(0, -1) & orthologs_85$sig_CyclinD1_MM %in% c(0, -1) & orthologs_85$sig_MIc_MM %in% c(0, -1) & orthologs_85$sig_MIc_MGUS %in% c(0, -1) & orthologs_85$sig_Mmset_MM %in% c(0, -1) & orthologs_85$sig_Mmset_MGUS %in% c(0, -1) & orthologs_85$sig_Trp53_MM %in% c(0, -1) & orthologs_85$sig_Trp53_MGUS %in% c(0, -1) & orthologs_85$sig_MGUS_HC %in% c(0, -1) & orthologs_85$sig_SMM_HC %in% c(0, -1) & orthologs_85$sig_MM_HC %in% c(0, -1),]
# remove rows with 0 for all human contrasts
orthologs_85_down <- orthologs_85_down[!(orthologs_85_down$sig_MGUS_HC == 0 & orthologs_85_down$sig_SMM_HC == 0 & orthologs_85_down$sig_MM_HC == 0),]


dim(orthologs_85_up) # 13
dim(orthologs_85_down) # 40


# Create gene sets
gene_lists <- list(
  orthologs_up = unique(orthologs_85_up$Mouse_ensembl_gene_id),
  orthologs_down = unique(orthologs_85_down$Mouse_ensembl_gene_id)
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
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
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
saveWorkbook(wb, file = "./ORA_orthologs.xlsx", overwrite = TRUE)

# Plotting
pdf("ORA_orthologs.pdf", width = 14, height = 7)
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


# CORCONDANT

# Other genes 

# Load overlaping orthologs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
orthologs_any <- read.table("orthologs_any.txt", sep = "\t", header = TRUE)

# ADD dynamics 
# Add differential results
orthologs_any <- merge(orthologs_any, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
orthologs_any <- merge(orthologs_any, human_dea_results[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
orthologs_any$sig_MM_HC <- human_dea_results_noscore[orthologs_any$Human_ensembl_gene_id,]$sig_MM_HC

length(unique(orthologs_any$Mouse_ensembl_gene_id)) # 341

# Keep rows only with 0 and 1
orthologs_any_up <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, 1) & orthologs_any$sig_CyclinD1_MM %in% c(0, 1) & orthologs_any$sig_MIc_MM %in% c(0, 1) & orthologs_any$sig_MIc_MGUS %in% c(0, 1) & orthologs_any$sig_Mmset_MM %in% c(0, 1) & orthologs_any$sig_Mmset_MGUS %in% c(0, 1) & orthologs_any$sig_Trp53_MM %in% c(0, 1) & orthologs_any$sig_Trp53_MGUS %in% c(0, 1) & orthologs_any$sig_MGUS_HC %in% c(0, 1) & orthologs_any$sig_SMM_HC %in% c(0, 1) & orthologs_any$sig_MM_HC %in% c(0, 1),]
#remove rows with 0 for all human contrasts
orthologs_any_up <- orthologs_any_up[!(orthologs_any_up$sig_MGUS_HC == 0 & orthologs_any_up$sig_SMM_HC == 0 & orthologs_any_up$sig_MM_HC == 0),]

# Keep rows only with 0 and -1
orthologs_any_down <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, -1) & orthologs_any$sig_CyclinD1_MM %in% c(0, -1) & orthologs_any$sig_MIc_MM %in% c(0, -1) & orthologs_any$sig_MIc_MGUS %in% c(0, -1) & orthologs_any$sig_Mmset_MM %in% c(0, -1) & orthologs_any$sig_Mmset_MGUS %in% c(0, -1) & orthologs_any$sig_Trp53_MM %in% c(0, -1) & orthologs_any$sig_Trp53_MGUS %in% c(0, -1) & orthologs_any$sig_MGUS_HC %in% c(0, -1) & orthologs_any$sig_SMM_HC %in% c(0, -1) & orthologs_any$sig_MM_HC %in% c(0, -1),]
# remove rows with 0 for all human contrasts
orthologs_any_down <- orthologs_any_down[!(orthologs_any_down$sig_MGUS_HC == 0 & orthologs_any_down$sig_SMM_HC == 0 & orthologs_any_down$sig_MM_HC == 0),]

length(unique(orthologs_any_up$Mouse_ensembl_gene_id)) # 31
length(unique(orthologs_any_down_heatmap.pdf$Mouse_ensembl_gene_id)) # 160


# Create gene sets
gene_lists <- list(
  orthologs_up = unique(orthologs_any_up$Mouse_ensembl_gene_id),
  orthologs_down = unique(orthologs_any_down$Mouse_ensembl_gene_id)
)


# Run ORA
ora_results <- list()
for (name in names(gene_lists)) {
  ensembl_ids <- gene_lists[[name]]
  ora_results[[name]] <- run_ora(ensembl_ids, name)
}


# Save Excel
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
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
saveWorkbook(wb, file = "./ORA_any_orthologs.xlsx", overwrite = TRUE)

# Plotting
pdf("ORA_any_orthologs.pdf", width = 14, height = 7)
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


# Human UP | Mouse DOWN
# Load overlaping orthologs
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
orthologs_any <- read.table("orthologs_any.txt", sep = "\t", header = TRUE)

# ADD dynamics 
# Add differential results
orthologs_any <- merge(orthologs_any, dea_results[,c("sig_CyclinD1_MGUS", "sig_CyclinD1_MM", "sig_MIc_MM", "sig_MIc_MGUS", "sig_Mmset_MM", "sig_Mmset_MGUS","sig_Trp53_MM", "sig_Trp53_MGUS")], by.x = "Mouse_ensembl_gene_id", by.y = "row.names")

# Add human results 
orthologs_any <- merge(orthologs_any, human_dea_results[,c("sig_MGUS_HC", "sig_SMM_HC")], by.x = "Human_ensembl_gene_id", by.y = "row.names")
orthologs_any$sig_MM_HC <- human_dea_results_noscore[orthologs_any$Human_ensembl_gene_id,]$sig_MM_HC

length(unique(orthologs_any$Mouse_ensembl_gene_id)) # 341

orthologs_up_down <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, -1) & orthologs_any$sig_CyclinD1_MM %in% c(0, -1) & orthologs_any$sig_MIc_MM %in% c(0, -1) & orthologs_any$sig_MIc_MGUS %in% c(0, -1) & orthologs_any$sig_Mmset_MM %in% c(0, -1) & orthologs_any$sig_Mmset_MGUS %in% c(0, -1) & orthologs_any$sig_Trp53_MM %in% c(0, -1) & orthologs_any$sig_Trp53_MGUS %in% c(0, -1) & orthologs_any$sig_MGUS_HC %in% c(0, 1) & orthologs_any$sig_SMM_HC %in% c(0, 1) & orthologs_any$sig_MM_HC %in% c(0, 1),]
#remove rows with 0 for all human contrasts
orthologs_up_down <- orthologs_up_down[!(orthologs_up_down$sig_MGUS_HC == 0 & orthologs_up_down$sig_SMM_HC == 0 & orthologs_up_down$sig_MM_HC == 0),]

length(unique(orthologs_up_down$Mouse_ensembl_gene_id)) # 56



orthologs_down_up <- orthologs_any[orthologs_any$sig_CyclinD1_MGUS %in% c(0, 1) & orthologs_any$sig_CyclinD1_MM %in% c(0, 1) & orthologs_any$sig_MIc_MM %in% c(0, 1) & orthologs_any$sig_MIc_MGUS %in% c(0, 1) & orthologs_any$sig_Mmset_MM %in% c(0, 1) & orthologs_any$sig_Mmset_MGUS %in% c(0, 1) & orthologs_any$sig_Trp53_MM %in% c(0, 1) & orthologs_any$sig_Trp53_MGUS %in% c(0, 1) & orthologs_any$sig_MGUS_HC %in% c(0, -1) & orthologs_any$sig_SMM_HC %in% c(0, -1) & orthologs_any$sig_MM_HC %in% c(0, -1),]
#remove rows with 0 for all human contrasts
orthologs_down_up <- orthologs_down_up[!(orthologs_down_up$sig_MGUS_HC == 0 & orthologs_down_up$sig_SMM_HC == 0 & orthologs_down_up$sig_MM_HC == 0),]

length(unique(orthologs_down_up$Mouse_ensembl_gene_id)) # 94



# Create gene sets
gene_lists <- list(
  orthologs_up_down = unique(orthologs_up_down$Mouse_ensembl_gene_id),
  orthologs_down_up = unique(orthologs_down_up$Mouse_ensembl_gene_id)
)


# Run ORA
ora_results <- list()
for (name in names(gene_lists)) {
  ensembl_ids <- gene_lists[[name]]
  ora_results[[name]] <- run_ora(ensembl_ids, name)
}


# Save Excel
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/DEA/")
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
saveWorkbook(wb, file = "./ORA_up_down_orthologs.xlsx", overwrite = TRUE)

# Plotting
pdf("ORA_up_down_orthologs.pdf", width = 14, height = 7)
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

