library(ComplexHeatmap)
library(dplyr)
library(readr)

# Load and prep data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/RNA/Count_matrix")
meta <- read.table("Col_Data.csv", sep = ",")
colnames(meta) <- meta[1,]
meta <- meta[-1,]
rownames(meta) <- meta$Sample_ID

# Load rna data for final samples 
load("norm_to_voom.RData")
meta <- meta[meta$Sample_ID %in% colnames(d1_norm$counts),]

# Change unknown sex into female
meta[meta$Sex == "Unknown", ]$Sex <- "female"



# Get ATAC data

# Load and prep data
setwd("/ibex/user/kurowsaa/Riney_project/Mouse/ATAC/Count_matrix")
met <- read.table("Col_Data.csv", sep = ",")
colnames(met) <- met[1,]
met <- met[-1,]
rownames(met) <- met$Sample_ID

# Load rna data for final samples 
load("norm_to_voom.RData")
met <- met[met$Sample_ID %in% colnames(d1_norm$counts),]

# Are all the ATAC within RNA ?
met$Sample_ID[!met$Sample_ID %in% meta$Sample_ID]

# Unify names 
meta$ATAC_ID <- NA
for (i in 1:nrow(met)) {
  # get ATAC ID
  atac <- met$Sample_ID[i]
  for (j in 1:nrow(meta)) {
    # get RNA ID
    rna <- meta$Sample_ID[j]
    # check if atac and rna are the same
    if (atac == rna) {
      # if they are the same, assign the ATAC ID to the RNA ID
      meta$ATAC_ID[j] <- atac
    }
  }
}

meta[meta$Stage == "Control",c("Sample_ID", "ATAC_ID")]
met[met$Stage == "Control", c("Sample_ID")]

meta["YFP_POOL_6_Control_Control_S27",]$ATAC_ID <- "Pool6_YFPCg1_Control_S6"
meta["YFP_POOL_7_Control_Control_S32",]$ATAC_ID <- "Pool7_YFPCg1_Control_S10"
meta["YFP_POOL_9_Control_Control_S35",]$ATAC_ID <- "Pool9_YFP_Control_S22"

meta[which(meta$Stage == "MGUS" & meta$Model == "CyclinD1_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MGUS" & met$Model == "CyclinD1_BIc"), c("Sample_ID")]

meta[which(meta$Stage == "MM" & meta$Model == "CyclinD1_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MM" & met$Model == "CyclinD1_BIc"), c("Sample_ID")]

meta["3069_CyclinD1BIc_MM_S18",]$ATAC_ID <- "3069_CyclinD1BIc_MM_S21"
meta["3162_CyclinD1BIc_MM_S17",]$ATAC_ID <- "3162_CyclinD1BIc_MM_S20"
meta["3661_CyclinD1BIc_MM_S22",]$ATAC_ID <- "3661_CyclinD1BIc_MM_S25"
meta["5278_CyclinD1BIc_MM_S19",]$ATAC_ID <- "5278_CyclinD1BIc_MM_S22"
meta["6575_CyclinD1BIc_MM_S20",]$ATAC_ID <- "6575_CyclinD1BIc_MM_S23"
meta["6844_CyclinD1BIc_MM_S21",]$ATAC_ID <- "6844_CyclinD1BIc_MM_S24"


meta[which(meta$Stage == "MGUS" & meta$Model == "MIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MGUS" & met$Model == "MIc"), c("Sample_ID")]

meta["8936_MIc_MGUS_S3",]$ATAC_ID <- "8936_MIC_MGUS_S4"
meta["8982_MIc_MGUS_S5",]$ATAC_ID <- "8982_MIC_MGUS_S3"
meta["8985_MIc_MGUS_S28",]$ATAC_ID <- "8985_MIC_MGUS_S1"
meta["8995_MIc_MGUS_S2",]$ATAC_ID <- "8995_MIC_MGUS_S2"
meta["SC80_MIc_MGUS_S7",]$ATAC_ID <- "SC80_MIC_MGUS_S27"

meta[which(meta$Stage == "MM" & meta$Model == "MIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MM" & met$Model == "MIc"), c("Sample_ID")]

meta["8476_MIc_Mieloma_S31",]$ATAC_ID <- "8476_MIC_Mieloma_S15"
meta["8582_MIc_MM_S17",]$ATAC_ID <- "8582_MIC_Mieloma_S5"
meta["8991_MIc_MM_S24",]$ATAC_ID <- "8991_MIC_MIeloma_S30"
meta["8994_MIc_MM_S9",]$ATAC_ID <- "8994_MIC_MIeloma_S28"
meta["SC30_MIc_MM_S21",]$ATAC_ID <- "SC30_MIC_Mieloma_S14"


meta[which(meta$Stage == "MGUS" & meta$Model == "Mmset_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MGUS" & met$Model == "Mmset_BIc"), c("Sample_ID")]

meta[which(meta$Stage == "MM" & meta$Model == "Mmset_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MM" & met$Model == "Mmset_BIc"), c("Sample_ID")]

meta["3086_MMSetB2IC_MM_S21",]$ATAC_ID <- "3086_MmsetBIc_MM_S17"

# Add from met  "3408_MmsetBIc_MM_S19" ,"3670_MmsetBIc_MM_S18",

meta[which(meta$Stage == "MGUS" & meta$Model == "Trp53_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MGUS" & met$Model == "Trp53_BIc"), c("Sample_ID")]

meta["8322_PBIc_MGUS_S20",]$ATAC_ID <- "8322_pB2IC_Cp_MGUS_S29"
meta["8326_PBIc_MGUS_S10",]$ATAC_ID <- "8326_pB2IC_Cp_MGUS_S32"
meta["8328_PBIc_MGUS_S29",]$ATAC_ID <- "8328_pB2IC_Cp_MGUS_S23"
meta["8330_PBIc_MGUS_S19",]$ATAC_ID <- "8330_pB2IC_Cp_MGUS_S24"
meta["8331_PBIc_MGUS_S12",]$ATAC_ID <- "8331_pB2IC_Cp_MGUS_S18"

meta[which(meta$Stage == "MM" & meta$Model == "Trp53_BIc"), c("Sample_ID", "ATAC_ID")]
met[c(met$Stage == "MM" & met$Model == "Trp53_BIc"), c("Sample_ID")]

meta["249_PBIc_MM_S4",]$ATAC_ID <- "249_pB2IC_Cp_Mieloma_S7"
meta["6174_PBIc_MM_S13",]$ATAC_ID <- "6174_pB2IC_Cp_Mieloma_S9"
meta["8323_PBIc_MM_S26",]$ATAC_ID <- "8323_pB2IC_Cp_Mieloma_S31"
meta["8327_PBIc_MM_S25",]$ATAC_ID <- "8327_pB2IC_Cp_Mieloma_S19"
meta["8329_PBIc_MM_S22",]$ATAC_ID <- "8329_pB2IC_Cp_Mieloma_S11"
meta["8332_PBIc_MM_S15",]$ATAC_ID <- "8332_pB2IC_Cp_Mieloma_S25"

# take first sample number before _ sign
# Are any of the samples in both omics?
table(gsub("_.*", "", meta[is.na(meta$ATAC_ID),]$Sample_ID)%in%gsub("_.*", "", met$Sample_ID))
table(gsub("_.*", "", meta[is.na(meta$ATAC_ID),]$Sample_ID)%in%gsub("_.*", "", met$Sample_ID))

toADD <- c("3408_MmsetBIc_MM_S19" ,"3670_MmsetBIc_MM_S18")
length(toADD) # 2

colnames(met)[14] <- "%PC BM "
met$ATAC_ID <- rownames(met)
colnames(meta)

met_atac <- met[,colnames(meta)]
dim(met_atac) # 49 x 14

colnames(met_atac) == colnames(meta)
dim(meta) # 72 x 14
rownames(met_atac) <- met_atac$ATAC_ID

table(meta$ATAC_ID %in% met_atac$ATAC_ID) # 47 TRUE
table(met_atac$ATAC_ID %in% meta$ATAC_ID) # 47 TRUE
met_atac[which(!met_atac$ATAC_ID %in% meta$ATAC_ID),]

# Add ATAC data to meta
metadata <- rbind(meta, met_atac[toADD,])
dim(metadata) # 74 x 14
# Save metadata of both omics
setwd("/ibex/user/kurowsaa/Riney_project/Mouse")
write.csv(metadata, "metadata.csv", row.names = FALSE)

# Select and reshape data
library(textshape)

# Create a separate row from each model from the Model row
metadata$CyclinD1_BIc <- ifelse(metadata$Model == "CyclinD1_BIc", "CyclinD1_BIc", NA)
metadata$MIc <- ifelse(metadata$Model == "MIc", "MIc", NA)
metadata$Mmset_BIc <- ifelse(metadata$Model == "Mmset_BIc", "Mmset_BIc", NA)
metadata$Trp53_BIc <- ifelse(metadata$Model == "Trp53_BIc", "Trp53_BIc", NA)

meta_toPlot <- metadata %>%
  column_to_rownames("Sample_ID") %>%
  t()

# Make OncoPrint
# Convert to list format expected by oncoPrint
meta_toPlot <- apply(meta_toPlot[c("CyclinD1_BIc", "MIc", "Mmset_BIc", "Trp53_BIc"),], c(1,2), function(x) ifelse(is.na(x), "", x))


# Top annotation: sex 
sex <- metadata$Sex
names(sex) <- metadata$Sample_ID

color_sex <- c("light blue", "pink", "aquamarine", "gray")
names(color_sex) <- c("male", "female", "mix", "Unknown")

top_anno <- HeatmapAnnotation(
  Sex = sex,
  annotation_name_side = "left",
  col = list(Sex = as.factor(color_sex))
)

# add stage annotation and devide samples per stage
stage <- metadata$Stage
color_stage <- c("#368716", "#D27D46", "#8D4E85") 
names(color_stage) <- c("Control", "MGUS", "MM")

bottom_anno <- HeatmapAnnotation(
  Stage = stage,
  annotation_name_side = "left",
  col = list(Stage = as.factor(color_stage))
)

setwd("/ibex/user/kurowsaa/Riney_project/Mouse")
# Plot
pdf("META_oncoprint.pdf", width = 12, height = 2.5)
oncoPrint(meta_toPlot,
  alter_fun = list(
    background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#F0F0F0", col = NA)),
    "CyclinD1_BIc" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#CAB2DC", col = NA)),
    "MIc" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#1F78B4", col = NA)),
    "Mmset_BIc" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#FB9A99", col = NA)),
    "Trp53_BIc" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#CAF2B0", col = NA))
  ),
  col = c("CyclinD1_BIc" = "#CAB2DC", 
          "MIc" = "#1F78B4",
          "Mmset_BIc" = "#FB9A99",
          "Trp53_BIc" = "#CAF2B0"),
  show_column_names = FALSE,
  # column title on the bottom
  row_labels = c("Cyclin D1", "c-MYC", "MMSET", "TP53"),
  column_title_side = "bottom",
  column_split = factor(metadata$Stage, levels = c("Control", "MGUS", "MM")),
  #column_title_gp = gpar(fontsize = 6),
  alter_fun_is_vectorized = TRUE,
  heatmap_legend_param = list(title = "Models", labels = c( "Cyclin D1 t(11;14)", "c-MYC", "MMSET t(4;14)", "TP53")),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  right_annotation = NULL
)
dev.off()

table(metadata$Stage)
# Control    MGUS      MM 
#       6      32      36 
table(metadata$Model)
    #  Control CyclinD1_BIc          MIc    Mmset_BIc    Trp53_BIc 
    #        6           18           14           19           17 