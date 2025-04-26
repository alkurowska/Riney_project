setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/03_chromVAR")
dea_results <- read.table("chromVAR_dea_results.txt", header = T, sep = "\t")
dea_res_noscore <- read.table("chromVAR_dea_results_noscore.txt", header = T, sep = "\t")

TFs <- as.data.frame(rownames(dea_results))
colnames(TFs) <- "JASPAR_ID"
TFs$TFs <- dea_results$TFs
rownames(TFs) <- TFs$JASPAR_ID

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
load("norm_to_voom.RData")
ids <- gene_anno[gene_anno$gene_name%in%TFs$TFs,]
TFs$gene_ID <- NA
for(i in 1:nrow(TFs)){
  TFs$gene_ID[i] <- ids[ids$gene_name == TFs$TFs[i],]$gene_id
}


TFs$chromVAR_logFC_MGUS_HC <- dea_results$logFC_MGUS_HC
TFs$chromVAR_p.adj_MGUS_HC <- dea_results$p.adj_MGUS_HC
TFs$chromVAR_logFC_SMM_HC <- dea_results$logFC_SMM_HC
TFs$chromVAR_p.adj_SMM_HC <- dea_results$p.adj_SMM_HC
TFs$chromVAR_logFC_MM_HC <- dea_res_noscore$logFC_MM_HC
TFs$chromVAR_p.adj_MM_HC <- dea_res_noscore$p.adj_MM_HC


# DEA
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
rna_results <- read.table("rna_dea_results.txt", header = T, sep = "\t")
TFs$DE_logFC_MGUS_HC <- rna_results[TFs$gene_ID,]$logFC_MGUS_HC
TFs$DE_p.adj_MGUS_HC <- rna_results[TFs$gene_ID,]$p.adj_MGUS_HC
TFs$DE_logFC_SMM_HC <- rna_results[TFs$gene_ID,]$logFC_SMM_HC
TFs$DE_p.adj_SMM_HC <- rna_results[TFs$gene_ID,]$p.adj_SMM_HC

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
rna_results <- read.table("rna_dea_results_noscore.txt", header = T, sep = "\t")
TFs$DE_logFC_MM_HC <- rna_results[TFs$gene_ID,]$logFC_MM_HC
TFs$DE_p.adj_MM_HC <- rna_results[TFs$gene_ID,]$p.adj_MM_HC

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/00_custom_pathways")
load("JASPAR2024.RData")


setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/02_OR")
MGUS_OR <- read.table("MGUS_HC_odds.txt", header = T, sep = "\t")
SMM_OR <- read.table("SMM_HC_odds.txt", header = T, sep = "\t")
MM_OR <- read.table("MM_HC_odds.txt", header = T, sep = "\t")

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/02_OR_inv")
MGUS_OR_inv <- read.table("MGUS_HC_odds.txt", header = T, sep = "\t")
SMM_OR_inv <- read.table("SMM_HC_odds.txt", header = T, sep = "\t")
MM_OR_inv <- read.table("MM_HC_odds.txt", header = T, sep = "\t")


TFs$OR_MGUS_HC <- MGUS_OR$odds_ratio
TFs$OR_p.adj_MGUS_HC <- MGUS_OR$p.adj
TFs$OR_SMM_HC <- SMM_OR$odds_ratio
TFs$OR_p.adj_SMM_HC <- SMM_OR$p.adj
TFs$OR_MM_HC <- MM_OR$odds_ratio
TFs$OR_p.adj_MM_HC <- MM_OR$p.adj
TFs$OR_inv_MGUS_HC <- MGUS_OR_inv$odds_ratio
TFs$OR_inv_p.adj_MGUS_HC <- MGUS_OR_inv$p.adj
TFs$OR_inv_SMM_HC <- SMM_OR_inv$odds_ratio
TFs$OR_inv_p.adj_SMM_HC <- SMM_OR_inv$p.adj
TFs$OR_inv_MM_HC <- MM_OR_inv$odds_ratio
TFs$OR_inv_p.adj_MM_HC <- MM_OR_inv$p.adj


# GSEA 
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/01_ES/MGUS_HC")
MGUS_abs <- read.table("MGUS_HC_GSEA_TF_abs.txt", header = T, sep = "\t")

TFs$abs_NES_MGUS_HC <- NA
TFs$abs_pval_MGUS_HC <- NA

TFs[MGUS_abs$pathway,]$abs_NES_MGUS_HC <- MGUS_abs$NES
TFs[MGUS_abs$pathway,]$abs_pval_MGUS_HC <- MGUS_abs$pval

MGUS_pos <- read.table("MGUS_HC_GSEA_TF_pos.txt", header = T, sep = "\t")

TFs$pos_NES_MGUS_HC <- NA
TFs$pos_pval_MGUS_HC <- NA

TFs[MGUS_pos$pathway,]$pos_NES_MGUS_HC <- MGUS_pos$NES
TFs[MGUS_pos$pathway,]$pos_pval_MGUS_HC <- MGUS_pos$pval


MGUS_neg <- read.table("MGUS_HC_GSEA_TF_neg.txt", header = T, sep = "\t")

TFs$neg_NES_MGUS_HC <- NA
TFs$neg_pval_MGUS_HC <- NA

TFs[MGUS_neg$pathway,]$neg_NES_MGUS_HC <- MGUS_neg$NES
TFs[MGUS_neg$pathway,]$neg_pval_MGUS_HC <- MGUS_neg$pval

setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/01_ES/SMM_HC")
SMM_abs <- read.table("SMM_HC_GSEA_TF_abs.txt", header = T, sep = "\t")

TFs$abs_NES_SMM_HC <- NA
TFs$abs_pval_SMM_HC <- NA

TFs[SMM_abs$pathway,]$abs_NES_SMM_HC <- SMM_abs$NES
TFs[SMM_abs$pathway,]$abs_pval_SMM_HC <- SMM_abs$pval

SMM_pos <- read.table("SMM_HC_GSEA_TF_pos.txt", header = T, sep = "\t")

TFs$pos_NES_SMM_HC <- NA
TFs$pos_pval_SMM_HC <- NA

TFs[SMM_pos$pathway,]$pos_NES_SMM_HC <- SMM_pos$NES
TFs[SMM_pos$pathway,]$pos_pval_SMM_HC <- SMM_pos$pval

SMM_neg <- read.table("SMM_HC_GSEA_TF_neg.txt", header = T, sep = "\t")

TFs$neg_NES_SMM_HC <- NA
TFs$neg_pval_SMM_HC <- NA

TFs[SMM_neg$pathway,]$neg_NES_SMM_HC <- SMM_neg$NES
TFs[SMM_neg$pathway,]$neg_pval_SMM_HC <- SMM_neg$pval


setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/01_ES/MM_HC")
MM_abs <- read.table("MM_HC_GSEA_TF_abs.txt", header = T, sep = "\t")

TFs$abs_NES_MM_HC <- NA
TFs$abs_pval_MM_HC <- NA

TFs[MM_abs$pathway,]$abs_NES_MM_HC <- MM_abs$NES
TFs[MM_abs$pathway,]$abs_pval_MM_HC <- MM_abs$pval

MM_pos <- read.table("MM_HC_GSEA_TF_pos.txt", header = T, sep = "\t")
TFs$pos_NES_MM_HC <- NA
TFs$pos_pval_MM_HC <- NA

TFs[MM_pos$pathway,]$pos_NES_MM_HC <- MM_pos$NES
TFs[MM_pos$pathway,]$pos_pval_MM_HC <- MM_pos$pval

MM_neg <- read.table("MM_HC_GSEA_TF_neg.txt", header = T, sep = "\t")
TFs$neg_NES_MM_HC <- NA
TFs$neg_pval_MM_HC <- NA
TFs[MM_neg$pathway,]$neg_NES_MM_HC <- MM_neg$NES
TFs[MM_neg$pathway,]$neg_pval_MM_HC <- MM_neg$pval


# Save
setwd("/ibex/user/kurowsaa/Riney_project/TF_analysis/04_ranking")
write.table(TFs, "TFs.txt", sep = "\t", col.names = T, row.names = F)