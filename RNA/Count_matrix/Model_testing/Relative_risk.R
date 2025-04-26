#####################
### RELATIVE RISK ###
#####################

library(reshape2)
library(scales)
library(ggplot2)
library(epitools)

# PLOTTING FUNCTION

RR_plot <- function(result, result.p.val, stage) {
  
  # Transform the data for plotting
  new_results <- setNames(melt(as.matrix(RR)), c('rows', 'columns', 'values'))
  new_results$rows <- factor(new_results$rows, levels=new_results$rows[order(unique(new_results$rows), decreasing = F)])
  p.values <- setNames(melt(as.matrix(RR_p.value)), c('rows', 'columns', 'values'))
  p.values$rows <- factor(p.values$rows, levels=p.values$rows[order(unique(p.values$rows), decreasing = F)])

  new_results$p.value <- p.values$values

  # PLOT
  p <- ggplot2::ggplot(new_results, aes(y=rows,x=columns,
                             fill = ifelse(p.value < 0.05, p.value, NA), 
                             label = round(values,4))) +
  geom_tile() + geom_text(size=4,colour = "black") +
  scale_fill_gradient2(low = "red", high = "blue",, midpoint = 0.025, name = "p value") +
  ylim(rev(levels(new_results$rows))) + 
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_discrete(expand = c(0, 0)) +
  facet_grid(~columns,scales="free") +
  labs(x = "", y= "") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=4),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 4),
        legend.title = element_text( size=4), legend.text=element_text(size=4),
        legend.box.margin = margin(0, 0, 0, 0)) + 
  theme_bw()

  ggsave(plot = p, dpi = 1800, filename = paste0(stage,"_RR.pdf"), width = 200, height = 100, units = "mm")
}

# Load the data
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
load("voom_to_dea.RData")

metadata <- fish_rna

setwd("/ibex/user/kurowsaa/Riney_project/RNA/Count_matrix/Model_testing/")

stage <- c("MM", "SMM", "MGUS")

for(i in 1:length(stage)){
  # Keep only specific stage samples
  meta_stage <- metadata[metadata$Stage == stage[i],]

  # Keep only relevant info 
  fish_stage <- meta_stage[,c("del17p", "del1p", "trans_4_14", "trans_14_16", "qall", "Sex")]

  # replace mepty values with "neutral" in all of the columns
  fish_stage <- as.data.frame(lapply(fish_stage, function(x) ifelse(x == "", "neutral", x)))
  colnames(fish_stage) <- c("17p del", "1p del", "t(4;14)", "t(14;16)", "1q aber", "Sex")

  # remove columns with all neutral values
  fish_stage <- fish_stage[,colSums(fish_stage != "neutral") > 0]

  # Calculate RR
  # Create an empty results data frame
  RR <- data.frame(matrix(NA, ncol=ncol(fish_stage), nrow=ncol(fish_stage)))
  colnames(RR) <- colnames(fish_stage)
  rownames(RR) <- colnames(fish_stage)

  # Create an empty results p.values data frame (fisher exact test)
  RR_p.value <- RR


  # RR
  for(j in 1:nrow(RR)){
    for(k in 1:ncol(RR)){
    RRtable <- table(fish_stage[,j],fish_stage[,k])
    res <- epitools::riskratio(RRtable,rev='both',method = 'wald')
    RR[j,k] <- res$measure[2,1]
    RR_p.value[j,k] <- res$p.value[2,1]
    }
  }
  
  #Save data
  write.table(RR, paste0(stage[i],"_RR.txt"), sep = '\t', col.names = T, row.names = T)
  write.table(RR_p.value, paste0(stage[i],"_RR_p.val.txt"), sep = '\t', col.names = T, row.names = T)

  RR_plot(RR, RR_p.value, stage[i])

}