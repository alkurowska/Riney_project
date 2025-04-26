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
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Human_validation/")
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

# Plot venndiagrams of intersection

library(ggVennDiagram)
library(ggplot2)
# Create a list of gene sets
genes_stes <- list(genes_all$Human_ensembl_gene_id, all_genes)

# Create a list of gene set names
genes_names <- c("Mouse homologs", "Human expressed")

# Create a Venn diagram
png("venn_diagram_all_genes.png", width=2200, height=1500, res=300)
ggVennDiagram(genes_stes, category.names = c("Mouse
homologs", "Human
expressed")) +
coord_flip()
dev.off()


genes_stes <- list(final$Human_ensembl_gene_id, final_genes)

# Create a list of gene set names
genes_names <- c("Mouse homologs", "Human expressed")

# Create a Venn diagram
png("venn_diagram_final_genes.png", width=2200, height=1500, res=300)
ggVennDiagram(genes_stes, category.names = c("Mouse
homologs", "Human
expressed")) +
coord_flip()
dev.off()