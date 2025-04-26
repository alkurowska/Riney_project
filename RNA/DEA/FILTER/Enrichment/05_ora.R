# Gene Ontology (GO) enrichment analysis 
# clusterProfiler package
# biological processes

# mydata - a dataframe with gene list and cluster assignment (up- or down- regualted)
# based on significantly differential genes, column names: geneID, cluster

ora <- function(mydata, contrast) {
    cluster_enrich_BP <- clusterProfiler::compareCluster(geneClusters = geneID ~ cluster, fun = "enrichGO", 
                                    OrgDb = org.Hs.eg.db, 
                                    data = mydata,
                                    ont = "BP", 
                                    universe = geneID,
                                    pvalueCutoff = 0.05,
                                    keyType = 'ENSEMBL')
    cluster_enrich_BP <- setReadable(cluster_enrich_BP, OrgDb = org.Hs.eg.db, keyType = "auto")
    write.table(cluster_enrich_BP,paste0(contrast ,"_cluster_enrichment_GO_BP.txt"), sep = "\t", col.names = T, row.names = F)
    saveRDS(cluster_enrich_BP, paste0(contrast ,"_cluster_enrichment_GO_BP.rds"))

    # plotting 
    BP <- dotplot(cluster_enrich_BP, font.size = 8, showCategory = 10, color = 'p.adjust', by='geneRatio', title = "Biological Processes")
    ggsave(plot = BP, width = 7, height = 8, dpi = 400, filename = paste0(contrast ,"_cluster_enrichment_GO_BP_dotplot_20.jpg"))
}



