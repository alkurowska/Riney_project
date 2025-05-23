###########################################
####     RNA-Seq pipeline - SINGLE     ####
####           04_upset_plots.R        ####
####              HUMAN                ####
###########################################

# LOAD Differential Results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER")
load("final_genes.RData")


# VARIANCE EXPLAINED 

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing")
load("varPart_data.RData")

library(variancePartition)
vp <- sortCols(varPart)
colnames(vp) <- c("GSVA", "Stage","t(14;16)",  "1q aber", "17p del", "1p del",  "Sex", "t(4;14)", "Residuals")

toPlot <- unique(c(MGUS_HC_filter, SMM_HC_filter, MM_HC_filter))
vp_dea <- vp[rownames(vp) %in% toPlot,]

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/FILTER/PLOTS")

png("var_frac_RNA.png", res = 300, width = 6, height = 4, units = "in")
plotVarPart(vp_dea)
dev.off()
