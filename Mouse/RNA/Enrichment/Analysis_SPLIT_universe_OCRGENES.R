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

# Create gene sets
gene_lists <- list(
CyclinD1_MGUS_up = CyclinD1_MGUS_HC_up,
CyclinD1_MGUS_down = CyclinD1_MGUS_HC_down,
Mmset_MGUS_up = Mmset_MGUS_HC_up,
Mmset_MGUS_down = Mmset_MGUS_HC_down,
Trp53_MGUS_up = Trp53_MGUS_HC_up,  
Trp53_MGUS_down = Trp53_MGUS_HC_down,
MIC_MGUS_up = MIc_MGUS_HC_up,
MIC_MGUS_down = MIc_MGUS_HC_down,
CyclinD1_MM_up = CyclinD1_MM_HC_up,
CyclinD1_MM_down = CyclinD1_MM_HC_down,
Mmset_MM_up = Mmset_MM_HC_up,
Mmset_MM_down = Mmset_MM_HC_down,
Trp53_MM_up = Trp53_MM_HC_up,
Trp53_MM_down = Trp53_MM_HC_down,
MIC_MM_up = MIc_MM_HC_up,
MIC_MM_down = MIc_MM_HC_down

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
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Enrichment/")
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


# UPset plot for up-regualted genes from ora_results list
# Load the UpSetR library
library(UpSetR)
up_names <- grep("_up$", names(ora_results), value = TRUE)

ora_results_up <- ora_results[up_names]

# Get pathways names 
CyclinD1_MM_up <-  ora_results_up$CyclinD1_MM_up$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
Mmset_MM_up <- ora_results_up$Mmset_MM_up$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
Trp53_MM_up <- ora_results_up$Trp53_MM_up$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
MIc_MM_up <- ora_results_up$MIC_MM_up$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
MIc_MGUS_up <- ora_results_up$MIC_MGUS_up$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)


sets <- unique(c(CyclinD1_MM_up, Mmset_MM_up, Trp53_MM_up, MIc_MM_up, MIc_MGUS_up))
# Create a binary matrix of membership

genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = genes,
  MIc_MGUS = genes %in% MIc_MGUS_up,
  CyclinD1_MM = genes %in% CyclinD1_MM_up,
  Mmset_MM = genes %in% Mmset_MM_up,
  Trp53_MM = genes %in% Trp53_MM_up,
  MIc_MM = genes %in% MIc_MM_up
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("MIc_MGUS", "CyclinD1_MM", "Mmset_MM", "Trp53_MM", "MIc_MM")

# Generate the UpSet plot

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Enrichment/")
png("upset_plot_up_pathways.png", width=3000, height=1500, res=300)
upset(upset_data, sets = c("MIc_MGUS", "CyclinD1_MM", "Mmset_MM", "Trp53_MM", "MIc_MM"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()


down_names <- grep("_down$", names(gene_lists), value = TRUE)

ora_results_down <- ora_results[down_names]

# Get pathways names for p.adj < 0.05
CyclinD1_MM_down <- ora_results_down$CyclinD1_MM_down$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
Mmset_MM_down <- ora_results_down$Mmset_MM_down$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
Trp53_MM_down <- ora_results_down$Trp53_MM_down$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
MIc_MM_down <- ora_results_down$MIC_MM_down$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)
MIc_MGUS_down <- ora_results_down$MIC_MGUS_down$GO_BP@result %>% filter(p.adjust < 0.05) %>% pull(Description)


ora_results_down$MIC_MGUS_down$GO_BP@result[6:10,]
sets <- unique(c(CyclinD1_MM_down, Mmset_MM_down, Trp53_MM_down, MIc_MM_down, MIc_MGUS_down))
# Create a binary matrix of membership

genes <- unique(unlist(sets))
membership_matrix <- data.frame(
  Gene = genes,
  MIc_MGUS = genes %in% MIc_MGUS_down,
  CyclinD1_MM = genes %in% CyclinD1_MM_down,
  Mmset_MM = genes %in% Mmset_MM_down,
  Trp53_MM = genes %in% Trp53_MM_down,
  MIc_MM = genes %in% MIc_MM_down
)

# Convert logical values to numeric for plotting
membership_matrix[, -1] <- lapply(membership_matrix[, -1], as.numeric)
# UpSetR requires the data in binary format without the gene column
upset_data <- membership_matrix[, -1]

colnames(upset_data) <- c("MIc_MGUS", "CyclinD1_MM", "Mmset_MM", "Trp53_MM", "MIc_MM")

# Generate the UpSet plot

setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Enrichment/")
png("upset_plot_down_pathways.png", width=3000, height=1500, res=300)
upset(upset_data, sets = c("MIc_MGUS", "CyclinD1_MM", "Mmset_MM", "Trp53_MM", "MIc_MM"),
      order.by = "freq", 
      main.bar.color = "steelblue")
dev.off()


up_data_frame <- matrix(NA, ncol=4, nrow=0)
colnames(up_data_frame) <- c("Category", "Term", "p.adjust", "GeneRatio")
for (group in names(ora_results_up)) {
  for (category in names(ora_results_up[[group]])) {
    result <- ora_results_up[[group]][[category]]
    if (is.null(result) || nrow(result) == 0) next
    df <- as.data.frame(result)
    df$Category <- group
    df$Term <- category
    up_data_frame <- rbind(up_data_frame, df)
  }
}

head(up_data_frame)
library(ggplot2)

# Step 1: Get top 10 per category
top_terms <- up_data_frame %>%
  group_by(Category) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

# Step 2: Find all extra rows in other categories where those terms also appear
extra_terms <- up_data_frame %>%
  filter(Description %in% top_terms$Description) %>%
  anti_join(top_terms, by = c("Category", "Description"))

# Step 3: Combine
toPlot <- bind_rows(top_terms, extra_terms)

# Step 4: Flag repeated terms
desc_cat_counts <- toPlot %>%
  distinct(Description, Category) %>%
  group_by(Description) %>%
  summarise(n_categories = n(), .groups = "drop") %>%
  mutate(repeated = n_categories > 1)

toPlot <- left_join(toPlot, desc_cat_counts, by = "Description")

# Prepare wide matrix of terms vs categories (fill NA with 0 or dummy)
mat <- toPlot %>%
  distinct(Category, Description) %>%
  mutate(val = 1) %>%
  tidyr::pivot_wider(names_from = Category, values_from = val, values_fill = 0) %>%
  tibble::column_to_rownames("Description") %>%
  as.matrix()

# Hierarchical clustering
clust <- hclust(dist(mat))
term_order <- clust$labels[clust$order]

# Add order to plot
toPlot$Description <- factor(toPlot$Description, levels = term_order)

# change column order 
toPlot$Category <- factor(toPlot$Category, levels = c("MIC_MGUS_up", "CyclinD1_MM_up", "Mmset_MM_up", "Trp53_MM_up", "MIC_MM_up"))

# example structure: ORA_df with columns Category, Term, p.adjust, GeneRatio
pdf("ORA_comapre_up.pdf", width = 10, height = 8)
ggplot(toPlot, aes(x = Category, y =  reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_text(size = 8)) +
  labs(size = "Gene Count", color = "p.adjust")
dev.off()



# Repeat for down 
down_names <- c("MIC_MGUS_down", "CyclinD1_MM_down", "Mmset_MM_down", "Trp53_MM_down", "MIC_MM_down")

ora_results_down <- ora_results[down_names]

down_data_frame <- matrix(NA, ncol=4, nrow=0)
colnames(down_data_frame) <- c("Category", "Term", "p.adjust", "GeneRatio")
for (group in names(ora_results_down)) {
  for (category in names(ora_results_down[[group]])) {
    result <- ora_results_down[[group]][[category]]
    if (is.null(result) || nrow(result) == 0) next
    df <- as.data.frame(result)
    df$Category <- group
    df$Term <- category
    down_data_frame <- rbind(down_data_frame, df)
  }
}

# Step 1: Get top 10 per category
top_terms <- down_data_frame %>%
  group_by(Category) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

# Step 2: Find all extra rows in other categories where those terms also appear
extra_terms <- up_data_frame %>%
  filter(Description %in% top_terms$Description) %>%
  anti_join(top_terms, by = c("Category", "Description"))

# Step 3: Combine
toPlot <- bind_rows(top_terms, extra_terms)

# Step 4: Flag repeated terms
desc_cat_counts <- toPlot %>%
  distinct(Description, Category) %>%
  group_by(Description) %>%
  summarise(n_categories = n(), .groups = "drop") %>%
  mutate(repeated = n_categories > 1)

toPlot <- left_join(toPlot, desc_cat_counts, by = "Description")

# Prepare wide matrix of terms vs categories (fill NA with 0 or dummy)
mat <- toPlot %>%
  distinct(Category, Description) %>%
  mutate(val = 1) %>%
  tidyr::pivot_wider(names_from = Category, values_from = val, values_fill = 0) %>%
  tibble::column_to_rownames("Description") %>%
  as.matrix()

# Hierarchical clustering
clust <- hclust(dist(mat))
term_order <- clust$labels[clust$order]

# Add order to plot
toPlot$Description <- factor(toPlot$Description, levels = term_order)

# change column order 
toPlot$Category <- factor(toPlot$Category, levels = c("MIC_MGUS_down", "CyclinD1_MM_down", "Mmset_MM_down", "Trp53_MM_down", "MIC_MM_down"))

# example structure: ORA_df with columns Category, Term, p.adjust, GeneRatio
pdf("ORA_comapre_down.pdf", width = 12, height = 8)
ggplot(toPlot, aes(x = Category, y =  reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text.y = element_text(size = 8)) +
  labs(size = "Gene Count", color = "p.adjust")
dev.off()