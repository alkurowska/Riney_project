# Gene Set Enrichment Analysis (GSEA)
# clusterProfiler package
# biological processes

# mydata - a named vector with logFC values from a given contrast ordered in a decreasing order
# names of the vector are gene IDs

gse <- function(mydata, contrast) {
    cluster_gse_BP <- clusterProfiler::gseGO(geneList = mydata, 
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        verbose = TRUE,
                        keyType = 'ENSEMBL')
    cluster_gse_BP <- setReadable(cluster_gse_BP, OrgDb = org.Hs.eg.db, keyType = "auto")
    write.table(cluster_gse_BP,paste0(contrast ,"_gse_GO_BP.txt"), sep = "\t")
    saveRDS(cluster_gse_BP, paste0(contrast ,"_gse_GO_BP.rds"))
        
    # plotting 

    # NES density score
    NES <- density(cluster_gse_BP$NES)
    png(paste0(contrast ,"_NES_density_plot.png"))
    plot(NES, main = "Density plot of NES")
    dev.off()

    BP <- dotplot(cluster_gse_BP, showCategory = 30, color = 'NES')
    ggsave(plot = BP, width = 10, height = 20, dpi = 400, filename = paste0(contrast ,"_gse_GO_BP_dotplot.jpg"))
}

