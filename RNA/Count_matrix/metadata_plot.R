library(ComplexHeatmap)
library(dplyr)
library(readr)

# Load and prep data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix")
meta <- read.table("MM_FISH_FINAL_R.txt", sep = "\t")
colnames(meta) <- meta[1,]
meta <- meta[-1,]

# remove NA
# Remove FISH data with NA values
meta <- meta[!is.na(meta$`1q amp`),] # 202

# Remove samples with empty RNA
meta <- meta[!(meta$`Sample RNAseq`==""),]
dim(meta) # 197

colnames(meta)[4] <- "17p del"
colnames(meta)[5] <- "1p del"
colnames(meta)[11] <- "1q aber"
colnames(meta)[6] <- "trans4_14"
colnames(meta)[7] <- "trans4_16"



alterations <- c("1q aber", "trans4_14", "trans4_16", "17p del", "1p del")

# Clean sample names
meta$`Sample RNAseq` <- make.unique(meta$`Sample RNAseq`)

# Select and reshape data
library(textshape)

mat <- meta %>%
  select(`Sample RNAseq`, all_of(alterations)) %>%
  mutate(across(-`Sample RNAseq`, ~ ifelse(.x != "neutral", .x, NA))) %>%
  column_to_rownames("Sample RNAseq") %>%
  t()

# Make OncoPrint
# Convert to list format expected by oncoPrint
mat <- apply(mat, c(1,2), function(x) ifelse(is.na(x), "", x))
# replace "t(4;14)" with "trans4_14"
mat <- apply(mat, c(1,2), function(x) ifelse(x == "t(4;14)" , "trans4_14", x))
# replace "t(4;16)" with "trans4_16"
mat <- apply(mat, c(1,2), function(x) ifelse(x == "t(14;16)", "trans14_16", x))



# Top annotation: sex and age
sex <- meta$Sex
names(sex) <- meta$`Sample RNAseq`

color_sex <- c("light blue", "pink")
names(color_sex) <- c("male", "female")

age <- as.numeric(meta$Age)
names(age) <- meta$`Sample RNAseq`

top_anno <- HeatmapAnnotation(
  Age = anno_barplot(age, border = FALSE, gp = gpar(fill = "gray")),
  Sex = sex,
  annotation_name_side = "left",
  col = list(Sex = as.factor(color_sex))
)

# add stage annotation and devide samples per stage
stage <- meta$Stage
color_stage <- c("#368716", "#D27D46","#D92F5E", "#8D4E85") 
names(color_stage) <- c("HC", "MGUS", "SMM", "MM")

bottom_anno <- HeatmapAnnotation(
  Stage = stage,
  annotation_name_side = "left",
  col = list(Stage = as.factor(color_stage))
)


# Plot
pdf("FISH_oncoprint.pdf", width = 12, height = 2.5)
oncoPrint(mat,
  alter_fun = list(
    background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#F0F0F0", col = NA)),
    "1q aber" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#C0AB52", col = NA)),
    "trans4_14" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#ED90A4", col = NA)),
    "trans14_16" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#C699E7", col = NA)),
    "17p del" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#4FBF85", col = NA)),
    "1p del" = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "#28BBD7", col = NA))
  ),
  col = c("1q aber" = "#C0AB52", "trans4_14" = "#ED90A4", 
          "trans14_16" = "#C699E7", "17p del" = "#4FBF85", "1p del" = "#28BBD7"),
  show_column_names = FALSE,
  # column title on the bottom
  row_labels = c("1q aberration", "t(4;14)", "t(14;16)", "del(17p)", "del(1p)"),
  column_title_side = "bottom",
  column_split = factor(meta$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
  #column_title_gp = gpar(fontsize = 6),
  alter_fun_is_vectorized = TRUE,
  heatmap_legend_param = list(title = "Alterations", labels = c( "1q aberration", "t(4;14)", "t(14;16)", "del(17p)", "del(1p)")),
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  right_annotation = NULL

)
dev.off()