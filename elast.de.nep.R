tig.20 <- subset(tig,idents="TIG1-20")
tig.50 <- subset(tig,idents="TIG1-50")

FeaturePlot(tig.20,features = c("ITGA8","DLGAP5"), reduction = "umap")+
  FeaturePlot(tig.50,features = c("ITGA8","DLGAP5"), reduction = "umap")


tig.20.de.markers <- FindMarkers(tig.20, ident.1 = "NEP", ident.2 = "CTL",group.by = "condition")
tig.50.de.markers <- FindMarkers(tig.50, ident.1 = "NEP", ident.2 = "CTL",group.by = "condition")

tig.de.markers <- FindMarkers(tig, ident.1 = "NEP", ident.2 = "CTL",group.by = "condition")
tig.de.markers <- tig.de.markers[tig.de.markers$p_val_adj<0.05 & abs(tig.de.markers$avg_log2FC)>1,]
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
#
library("org.Hs.eg.db")
library(clusterProfiler)

tig.de.markers <- ms_ref[ms_ref$gene_short_name %in% rownames(tig.de.markers),]

ego_result <- enrichGO(gene          = tig.de.markers$entrezgene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))
