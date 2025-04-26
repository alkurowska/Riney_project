# Setting working directory
set.seed(123)
setwd("/ibex/user/kurowsaa/RINEY/human/Stratification/")

####################### LOAD THE DATA #######################

# Gene annotation 
gene_anno <- readRDS("Complete_GTF.RDS") #60666 x 12

### STEP 2: REMOVE Ig genes with chain IGH, IGK and IGL 
head(gene_anno)
Igs <- gene_anno[grep("^IG", gene_anno$gene_biotype),]
gene_annos <- gene_anno[!gene_anno$gene_id %in% Igs$gene_id,] #60252    12 #Jchain
dim(gene_annos)

### STEP 3: Are there any other genes?
levels(as.factor(gene_annos$gene_biotype))

Igs_other <- gene_annos[grep("^IGL|^IGH|^IGK", gene_annos$gene_name),]
gene_annos <- gene_annos[!gene_annos$gene_id %in% Igs_other$gene_id,] #60243    12
dim(gene_annos) 


rest <- gene_annos[grep("^IG", gene_annos$gene_name),] # 53
sort(unique(rest$gene_name))
# IGBP - immunoglobulin binding protein - REMOVE
# IGDCC - immunoglobulin superfamily, DCC subclass - REMOVE
# IGF1 - insulin like growth factor 1 - KEEP
# IGFN1 - REMOVE
# IGIP - IgA inducing protein - REMOVE
# IGSF - immunoglobulin super family - REMOVE

toRemove <- gene_annos[grep("^IGBP|^IGDCC|^IGFN1|^IGIP|^IGSF", gene_annos$gene_name),]
gene_annot <- gene_annos[!gene_annos$gene_id %in% toRemove$gene_id,] #60217    12
dim(gene_annot)

gene_annot[grep("^IG", gene_annot$gene_name),] # 53


### STEP 4: Save the gene annotation
setwd("/ibex/user/kurowsaa/Riney_project/RNA/QC")
saveRDS(gene_annot, "Gene_anno_noIGs.RDS")