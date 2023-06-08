
# Filtering
Asso_LFD <- rownames(Asso_LFD_R1[ Asso_LFD_R1$FDR <= 0.05 & rownames(Asso_LFD_R1) %in% rownames(Asso_LFD_R2[Asso_LFD_R2$FDR <= 0.05,]),])
Asso_HFD <- rownames(Asso_HFD_R1[ Asso_HFD_R1$FDR <= 0.05 & rownames(Asso_HFD_R1) %in% rownames(Asso_HFD_R2[Asso_HFD_R2$FDR <= 0.05,]),])

# Combine
Ass <- unique(sort(c(Asso_LFD, Asso_HFD)))

## Heatmap of genes associated with pseudo-time in LFD ------------------------------------------
## Predict smoothend expression 
Smooth_R1 <- predictSmooth(LFD_sce_R1, gene = Asso_LFD, tidy=F, n=100)
Smooth_R2 <- predictSmooth(LFD_sce_R2, gene = Asso_LFD, tidy=F, n=100)

# Average across replicates and scale
Smooth <- Smooth_R1
for (i in 1:nrow(Smooth)) { Smooth[i,] <- colMeans(rbind(Smooth_R1[i,], Smooth_R2[i,])) }
Smooth <- t(scale(t(Smooth)))

# Seriate the results
Smooth <- Smooth[ get_order(seriate(Smooth, method="PCA_angle")),]

# Create heatmaps (in 4 visually defined groups)
HTM1 <- Heatmap(Smooth[1:210,], cluster_columns=F, cluster_rows=F)
HTM2 <- Heatmap(Smooth[211:320,], cluster_columns=F, cluster_rows=F)
HTM3 <- Heatmap(Smooth[321:800,], cluster_columns=F, cluster_rows=F)
HTM4 <- Heatmap(Smooth[801:1255,], cluster_columns=F, cluster_rows=F)

# Plot the heatmap
HTM3 %v% HTM4 %v% HTM1 %v% HTM2

### Pathway analysis -----------------------------------------------
# Extract genes for each group
C1 <- rownames(Smooth[1:210,])
C2 <- rownames(Smooth[211:320,])
C3 <- rownames(Smooth[321:800,])
C4 <- rownames(Smooth[801:1255,])

# Setup
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Mus musculus", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) # TERM2NAME

# Convert gene names
C1_Entrez <- clusterProfiler::bitr(C1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C2_Entrez <- clusterProfiler::bitr(C2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C3_Entrez <- clusterProfiler::bitr(C3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
C4_Entrez <- clusterProfiler::bitr(C4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Do the pathway analysis
wiki_clust1 <- clusterProfiler::enricher(C1_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1,  TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust2 <- clusterProfiler::enricher(C2_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust3 <- clusterProfiler::enricher(C3_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
wiki_clust4 <- clusterProfiler::enricher(C4_Entrez[[2]], pAdjustMethod = "fdr", pvalueCutoff = 1,qvalueCutoff = 1, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

# Adding gene symbols to the resulting pathway file
wiki_clust1 <- as.data.frame(DOSE::setReadable(wiki_clust1, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust2 <- as.data.frame(DOSE::setReadable(wiki_clust2, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust3 <- as.data.frame(DOSE::setReadable(wiki_clust3, org.Mm.eg.db, keyType = "ENTREZID"))
wiki_clust4 <- as.data.frame(DOSE::setReadable(wiki_clust4, org.Mm.eg.db, keyType = "ENTREZID"))

## Barplot of top 3 pathways
# Cluster 1: Fasn, Scd1, Acaca, Pparg, Fabp4, Lep, Lipe, Plin1, Cd36
pdf("Wikipathway_Cluster1.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust1[1:3,"pvalue"]), names=wiki_clust1[1:3,"Description"], las=2)
dev.off()

# Cluster 2: Foxo1, Pten, Igf1, Akt2, Gsk3b, Ehd2
pdf("Wikipathway_Cluster2.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust2[1:3,"pvalue"]), names=wiki_clust2[1:3,"Description"], las=2)
dev.off()

# Cluster 3: Atp5a1, Rpl37, Rps2, Cov5a, Ndufa4,  Ubcrb
pdf("Wikipathway_Cluster3.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust3[1:3,"pvalue"]), names=wiki_clust3[1:3,"Description"], las=2)
dev.off()

# Cluster 4: Egfr, Pdgfrb, Pdgfra, Eps8, Eps15, Spry2
pdf("Wikipathway_Cluster4.pdf", width=10, height=10, useDingbats=F)
barplot(-log10(wiki_clust1[1:3,"pvalue"]), names=wiki_clust1[1:3,"Description"], las=2)
dev.off()