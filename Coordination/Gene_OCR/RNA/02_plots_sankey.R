####################
### sankey plots ###
####################

#### LOAD THE DATA ####
# Load the DEA results

# Load the DEA results
setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/GSVA")
dea_results <- read.table("rna_dea_results.txt", header=TRUE, row.names=1, sep="\t")

setwd("/ibex/user/kurowsaa/Riney_project/RNA/DEA/No_score")
dea_res_noscore <- read.table("rna_dea_results_noscore.txt", header=TRUE, row.names=1, sep="\t")

# Load significant genes 
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR")
load("coordination_genes.RData")

genes_toPlot <- unique(c(MGUS_HC, SMM_HC, MM_HC))

toPlot <- dea_results[genes_toPlot, c("sig_MGUS_HC","sig_SMM_HC")]
toPlot <- cbind(toPlot, dea_res_noscore[genes_toPlot, c("sig_MM_HC")])

#### PREP THE DATA ####
colnames(toPlot) <- c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC")

# Change 0s to "non differential"
toPlot[toPlot == 0] <- "Non-differential"
#Change -1 to "Down"
toPlot[toPlot == -1] <- "Down-regulated"
#Change 1 to "Up"
toPlot[toPlot == 1] <- "Up-regulated"

toPlot$Gene <- rownames(toPlot)

# Sanky plot
library(ggplot2)
library(ggalluvial)

# Reshape for ggalluvial
genes_long <- tidyr::pivot_longer(
  toPlot,
  cols = c(MGUS_vs_HC, SMM_vs_HC, MM_vs_HC),
  names_to = "Condition",
  values_to = "Category"
)

# Order conditions
genes_long$Condition <- factor(genes_long$Condition, 
                                levels = c("MGUS_vs_HC", "SMM_vs_HC", "MM_vs_HC"))

# Plot
setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")
png("rna_dea_sanky.png", width=1500, height=1000, res=300)
ggplot(genes_long, aes(
  x = Condition, stratum = Category, alluvium = Gene,
  fill = Category, label = Category
)) +
  geom_flow(stat = "alluvium") +
  geom_stratum() +
  theme_classic() +
  labs(x = NULL, y = "Gene Count") +
  scale_fill_manual(values = c("Up-regulated" = "#D7342A", "Down-regulated" = "#4475B3", "Non-differential" = "#F0F0F0"))
dev.off()


# # Plot a sankey plot by categorizing genes based on overlap with the contrasts

# library(networkD3)
# # Define node names (unique categories across conditions)
# nodes <- data.frame(name = c(
#   "MGUS_Down-regulated", "MGUS_Up-regulated", "MGUS_Non-differential",
#   "SMM_Down-regulated", "SMM_Up-regulated", "SMM_Non-differential",
#   "MM_Down-regulated", "MM_Up-regulated", "MM_Non-differential"
# ))

# # Convert categories to numeric node indices
# get_node_index <- function(category, condition) {
#   paste(condition, category, sep = "_") %>%
#     match(nodes$name) - 1
# }

# # Convert to long format (for links)
# library(dplyr)

# toPlot_long <- toPlot %>%
#   mutate(
#     MGUS_category = paste0("MGUS_", MGUS_vs_HC),
#     SMM_category = paste0("SMM_", SMM_vs_HC),
#     MM_category = paste0("MM_", MM_vs_HC)
#   )

# # Create links from MGUS → SMM, SMM → MM
# links <- toPlot_long %>%
#   group_by(MGUS_category, SMM_category, MM_category) %>%
#   summarise(weight = n(), .groups = "drop") %>%
#   mutate(
#     source1 = get_node_index(sub("MGUS_", "", MGUS_category), "MGUS"),
#     target1 = get_node_index(sub("SMM_", "", SMM_category), "SMM"),
#     source2 = get_node_index(sub("SMM_", "", SMM_category), "SMM"),
#     target2 = get_node_index(sub("MM_", "", MM_category), "MM")
#   ) %>%
#   select(source1, target1, source2, target2, weight)

# # Convert to networkD3 format (first set of transitions)
# links_sankey <- links %>%
#   select(source1, target1, weight) %>%
#   rename(source = source1, target = target1) %>%
#   bind_rows(links %>%
#               select(source2, target2, weight) %>%
#               rename(source = source2, target = target2))


# # Assign color to links based on "Special" column
# color_map <- JS("d3.scaleOrdinal()
#                   .domain(['MGUS_Up-regulated', 'MGUS_Non-differential', 'MGUS_Down-regulated', 
#                            'SMM_Up-regulated', 'SMM_Non-differential', 'SMM_Down-regulated',
#                            'MM_Up-regulated', 'MM_Non-differential', 'MM_Down-regulated'])
#                   .range(['#D7342A', '#F0F0F0', '#4475B3',
#                           '#D7342A', '#F0F0F0', '#4475B3',
#                           '#D7342A', '#F0F0F0', '#4475B3'])")

# library(htmlwidgets)

# # Save the interactive Sankey diagram as an HTML file
# s <- sankeyNetwork(
#   Links = links_sankey,
#   Nodes = nodes,
#   Source = "source",
#   Target = "target",
#   Value = "weight",
#   NodeID = "name",
#   #LinkGroup = "transition_type",  # This applies the colors to links!
#   colourScale = color_map,  # Define the link colors
#   sinksRight = TRUE,
#   fontSize = 0,
#   nodeWidth = 30
# )

# require(htmlwidgets)
# library(webshot)

# saveWidget(s, file="rna_dea_sankey_network.html")
# # webshot("rna_dea_sankey_network.html" , "rna_dea_sankey_network.pdf", delay = 0.2)

# # # Make a webshot in png : Low quality - but you can choose shape
# # webshot("rna_dea_sankey_network.html" , "output.png", delay = 0.2 , cliprect = c(440, 0, 1000, 10))

# # Assign colors to links based on transitions




# # Assign transistions
# links_sankey <- links_sankey %>%
#   mutate(
#     transition_type = case_when(
#       (source == 1 & target == 4) | (source == 4 & target == 7) ~ "Stable_Up",
#       (source == 0 & target == 3) | (source == 3 & target == 6) ~ "Stable_Down",
#       (source == 2 & target == 5) | (source == 5 & target == 8) ~ "Stable_Non",
#       TRUE ~ "Other"
#     )
#   )

# colors_map <- JS("d3.scaleOrdinal()
#                   .domain(['MGUS_Up-regulated', 'MGUS_Non-differential', 'MGUS_Down-regulated', 
#                            'SMM_Up-regulated', 'SMM_Non-differential', 'SMM_Down-regulated',
#                            'MM_Up-regulated', 'MM_Non-differential', 'MM_Down-regulated',
#                            'Stable_Up', 'Stable_Down', 'Stable_Non', 'Other'])
#                   .range(['#D7342A', '#F0F0F0', '#4475B3',
#                           '#D7342A', '#F0F0F0', '#4475B3',
#                           '#D7342A', '#F0F0F0', '#4475B3',
#                           '#CC6666', '#6699CC', '#F0F0F0', '#999999'])")

# ss <- sankeyNetwork(
#   Links = links_sankey,
#   Nodes = nodes,
#   Source = "source",
#   Target = "target",
#   Value = "weight",
#   NodeID = "name",
#   NodeGroup = "name",
#   LinkGroup = "transition_type",  # Apply colors based on event classification
#   colourScale = colors_map,  
#   sinksRight = TRUE,
#   fontSize = 0,
#   nodeWidth = 30
# )

# saveWidget(ss, file="rna_dea_sankey_network_transition.html")

# library(dplyr)
# library(tidyr)

# toPlot_long$MGUS_SMM_trans <- "transition"
# toPlot_long$SMM_MM_trans <- "transition"

# toPlot_long$MGUS_SMM <- ifelse(gsub("MGUS_","", toPlot_long$MGUS_category) == gsub("SMM_","", toPlot_long$SMM_category), 
#                                paste0("Stable_",gsub("MGUS_","",toPlot_long$MGUS_category)),
#                                "Other")

# toPlot_long$SMM_MM <- ifelse(gsub("SMM_","", toPlot_long$SMM_category) == gsub("MM_","", toPlot_long$MM_category), 
#                                paste0("Stable_",gsub("SMM_","",toPlot_long$SMM_category)),
#                                "Other")

# # Plot % barplot based on the transitions toPlot_long$SMM_MM or toPlot_long$MGUS_SMM
# setwd("/ibex/user/kurowsaa/Riney_project/Coordination/Gene_OCR/RNA/PLOTS")

# category_counts <- toPlot_long %>%
#   group_by(MGUS_SMM) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   mutate(percent = n / sum(n) * 100)

# colnames(category_counts) <- c("Category", "n", "percent")
# category_counts <- as.data.frame(category_counts)
# category_counts$Category <- factor(category_counts$Category, levels = c("Stable_Down-regulated", "Stable_Up-regulated", "Stable_Non-differential", "Other"))


# png("rna_dea_category_barplot_MGUS_SMM.png", width=1500, height=100, res=300)
# ggplot(category_counts, aes(x = 1, y = percent, fill = Category)) +
#   geom_bar(stat = "identity", width = 1) +
#   # Add % values 
#   geom_text(aes(label = paste0(round(percent,1), "%")), 
#             position = position_stack(vjust = 0.5), 
#             size = 2, color = "black") +  # Small font, white text inside bars
#   # keep the order as levels of the factor
#   scale_fill_manual(values = c( "Stable_Down-regulated" = "#6699CC", "Stable_Non-differential" = "#F0F0F0", "Stable_Up-regulated" = "#CC6666", "Other" = "#999999")) +
#   # flip -90 degrees
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = "none", panel.background = element_blank())
# dev.off()



# category_counts <- toPlot_long %>%
#   group_by(SMM_MM) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   mutate(percent = n / sum(n) * 100)

# colnames(category_counts) <- c("Category", "n", "percent")
# category_counts <- as.data.frame(category_counts)
# category_counts$Category <- factor(category_counts$Category, levels = c("Stable_Down-regulated", "Stable_Up-regulated", "Stable_Non-differential", "Other"))

# png("rna_dea_category_barplot_SMM_MM.png", width=1500, height=100, res=300)
# ggplot(category_counts, aes(x = 1, y = percent, fill = Category)) +
#   geom_bar(stat = "identity", width = 1) +
#   # Add % values 
#   geom_text(aes(label = paste0(round(percent,1), "%")), 
#             position = position_stack(vjust = 0.5), 
#             size = 2, color = "black") +  # Small font, white text inside bars
#   # keep the order as levels of the factor
#   scale_fill_manual(values = c( "Stable_Down-regulated" = "#6699CC", "Stable_Non-differential" = "#F0F0F0", "Stable_Up-regulated" = "#CC6666", "Other" = "#999999")) +
#   # flip -90 degrees
#   coord_flip() +
#   theme_void() +
#   theme(legend.position = "none", panel.background = element_blank())
# dev.off()

