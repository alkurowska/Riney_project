# Setting working directory
set.seed(123)

#Libraries
library(viridis)
library(ggplot2)
library(reshape2)
library(scales)

########################### PLOTTING ##########################


# Create functions for plotting 

# Function
anov_plot <- function(result, stage) {
  # Prepare data for ploting
  new_results <- setNames(melt(as.matrix(results)), c('PC', 'variable', 'values'))
  new_results$PC <- factor(new_results$PC, levels=new_results$PC[order(unique(new_results$PC), decreasing = F)])
  new_results$values <- as.numeric(new_results$values)

    p <- ggplot2::ggplot(new_results, aes(y=PC,x=variable,
                             fill = ifelse(values < 0.05, values, NA), label = round(values,4))) +
    geom_tile() + geom_text(size=2,colour = "black") +
    scale_fill_gradient2(low = "red", high = "blue",, midpoint = 0.025, name = "p value") +
    ylim(rev(levels(new_results$PC))) + 
    scale_x_discrete(expand = c(0, 0)) +
    #scale_y_discrete(expand = c(0, 0)) +
    facet_grid(~variable,scales="free") +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 5),
        legend.title = element_text( size=4), legend.text=element_text(size=4),
        legend.box.margin = margin(0, 0, 0, 0))
    
    ggsave(plot = p, dpi = 1800, filename = paste0(stage, "_anova.pdf"), width = 125, height = 100, units = "mm")
}




###########################  ANNOVA  ##########################

####################### DATA PREPARATION #######################

### STEP 1: LOAD DATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix")
log_cpm <- read.table("log_cpm.txt", header = T, row.names = 1)
colnames(log_cpm) <- gsub("X", "", colnames(log_cpm))

### STEP 2: LOAD METADATA
setwd("/ibex/user/kurowsaa/Riney_project/ATAC/DEA/GSVA")
load("voom_to_dea.RData")

setwd("/ibex/user/kurowsaa/Riney_project/ATAC/Count_matrix/Model_testing/Anova")

stage <- c("MGUS", "SMM", "MM")

for(i in 1:length(stage)){
  ## Perform PCA only for samples in one stage
  meta_stage <- fish_atac[fish_atac$Stage == stage[i],]
  norm_stage <- norm_data[,meta_stage$`Sample ATACseq`]
  pca_rna <- prcomp(t(norm_stage))

  # Plot PCs variance
  var_prcomp <- pca_rna$sdev^2
  pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), pc=c(1:length(var_prcomp)))
  pcvar$Y <- "ATACseq"
  pcvar$X <- paste0("PC", 1:nrow(pcvar))

  toPlot <- pcvar[1:10,]
  toPlot$X <- factor(toPlot$X, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
  toPlot$var <- toPlot$var*100

  p <- ggplot(toPlot, aes(X, Y, fill= var)) + 
    geom_tile() + 
    scale_fill_viridis(discrete=FALSE) +
    labs(x = "", y= "") + 
    coord_fixed(ratio=1) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))

  ggsave(plot = p, dpi = 1800, filename = paste0(stage[i],"_pcvar.pdf"), width = 125, height = 100, units = "mm")


  # ANOVA
  # pca_results 
  pc_scores  <- as.data.frame(pca_rna$x)

  # add metadata variable to the dataframe
  fish_stage <- meta_stage[,c("del17p", "del1p", "trans_4_14", "trans_14_16", "qall", "Sex", "gsva")]
  colnames(fish_stage) <- c("17p del", "1p del", "t(4;14)", "t(14;16)", "1q aber", "Sex", "GSVA")

  # Replace every empty cell with "neutral"
  fish_stage <- as.data.frame(lapply(fish_stage, function(x) ifelse(x == "", "neutral", x)))

  # check if all the values in the metadata are "neutral" - remove the columns with only neutral values
  fish_stage <- fish_stage[, colSums(fish_stage == "neutral") != nrow(fish_stage)]

  df <- cbind(pc_scores, fish_stage)

  # specify nr of components 
  num_pcs <- 10

  #results dataframe
  results <- data.frame(matrix(NA, ncol=length(fish_stage), nrow=10))

  metadata_var <- colnames(fish_stage)
  colnames(results) <- metadata_var

  colnames(df) <- c(colnames(pc_scores), colnames(fish_stage))

  for(j in 1:num_pcs){
    #Choose principal component
    pc_name <- paste0("PC", j)
    rownames(results)[j] <- pc_name
  
    # perform anova for each PC with metadata_var
    for(k in 1:length(metadata_var)){
      variable <- metadata_var[k]
        # Check if the variable is numeric
        if (is.numeric(df[[variable]])) {
          model <- lm(as.formula(paste0(pc_name, " ~ ", variable)), data = df)
          anova_result <- model
          # print summary 
          cat("Summary for", pc_name, ":\n")
          print(summary(anova_result))
    
          #Extract and print p-value
          p_value <- summary(anova_result)$coefficients[2, 4]
          cat("P-value for", pc_name, ":", p_value, "\n\n")

        } else {
          model <- aov(get(pc_name) ~ df[,variable], data = df)
          anova_result <- model
          # print summary 
          cat("Summary for", pc_name, ":\n")
          print(summary(anova_result))
    
          #Extract and print p-value
          p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
          cat("P-value for", pc_name, ":", p_value, "\n\n")
        }    

      results[pc_name,variable] <- p_value
    }
  }
  
  # remove X from the begning of column names
  colnames(results) <- gsub("X", "", colnames(results))
  colnames(results) <- gsub("17p.del", "17p del", colnames(results)) 
  colnames(results) <- gsub("1p.del", "1p del", colnames(results))
  colnames(results) <- gsub("t.4.14.", "t(4;14)", colnames(results))
  colnames(results) <- gsub("t.14.16.", "t(14;16)", colnames(results))
  colnames(results) <- gsub("1q.aber", "1q aber", colnames(results))
  # save results
  write.table(results, paste0(stage[i], "_anova_results.txt"))

  # Plot ANOVA results
  anov_plot(results, stage[i])

}



########################### ADJUSTED P-VALUES ##########################
# FUNCTION 
#Adjust for multiple testing for each column
anov_plot_adj <- function(results, stage) {
    results_p.adj <- results

    for(l in 1:nrow(results)){
       results_p.adj[l,] <- p.adjust(results[l,], method = "fdr")
       write.table(results_p.adj, paste0(stage, "_anova_results_p.adj.txt"))
    }

    new_results <- setNames(melt(as.matrix(results_p.adj)), c('PC', 'variable', 'values'))
    new_results$PC <- factor(new_results$PC, levels=new_results$PC[order(unique(new_results$PC), decreasing = F)])
    new_results$values <- as.numeric(new_results$values)

    p <- ggplot(new_results, aes(y=PC,x=variable, fill = ifelse(values < 0.05, values, NA), label = round(values,4))) +
      geom_tile() + geom_text(size=2,colour = "black") +
      scale_fill_gradient2(low = "red", high = "blue",, midpoint = 0.025, name = "p.adj") +
      ylim(rev(levels(new_results$PC))) + 
      scale_x_discrete(expand = c(0, 0)) +
      #scale_y_discrete(expand = c(0, 0)) +
      facet_grid(~variable,scales="free") +
      theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 5),
        legend.title = element_text( size=4), legend.text=element_text(size=4),
        legend.box.margin = margin(0, 0, 0, 0))

    ggsave(plot = p, dpi = 1800, filename = paste0(stage, "_anova_p.adj.pdf"), width = 125, height = 100, units = "mm")
  }


stage = c("MGUS","SMM", "MM")
# RUN 
  for(i in 1:length(stage)){
    # Load results
    results <- read.table(paste0(stage[i], "_anova_results.txt"), header = T, row.names = 1)
    # remove X from the begning of column names
    colnames(results) <- gsub("X", "", colnames(results))
    colnames(results) <- gsub("17p.del", "17p del", colnames(results)) 
    colnames(results) <- gsub("1p.del", "1p del", colnames(results))
    colnames(results) <- gsub("t.4.14.", "t(4;14)", colnames(results))
    colnames(results) <- gsub("t.14.16.", "t(14;16)", colnames(results))
    colnames(results) <- gsub("1q.aber", "1q aber", colnames(results))
    anov_plot_adj(results, stage[i])
  }

