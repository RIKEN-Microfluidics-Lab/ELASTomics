library(ggplot2)
library(ggrepel)
library(Seurat)
library(pheatmap)
library(destiny)
#https://satijalab.org/seurat/articles/integration_introduction.html
tig.list <- SplitObject(tig, split.by = "condition")
tig.list <- lapply(X = tig.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = tig.list)

tig.anchors <- FindIntegrationAnchors(object.list = tig.list, anchor.features = features)
tig.combined <- IntegrateData(anchorset = tig.anchors)
DefaultAssay(tig.combined) <- "integrated"
tig.combined <- ScaleData(tig.combined, verbose = FALSE)
tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindClusters(tig.combined, resolution = 0.1)
tig.combined <- RenameIdents(tig.combined, `0` = "TIG1-50", `1` = "TIG1-20")
tig.combined[["age"]] <- Idents(tig.combined)

DimPlot(tig, reduction = "pca",group.by = "condition")+
  DimPlot(tig, reduction = "pca")+
  DimPlot(tig.combined, reduction = "pca", group.by = "condition")+
  DimPlot(tig.combined, reduction = "pca", label = TRUE, repel = TRUE)

FeaturePlot(tig.combined,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","dtd_FLD040"), reduction = "umap")



# remove low quality samples
# tig.combined <- FindClusters(tig.combined, resolution = 0.2)
# DimPlot(tig.combined)
# FeatureScatter(tig.combined,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")


tig.combined <- NormalizeData(tig.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig.combined.nep <- subset(tig.combined,subset=condition=="NEP")
tig.combined.ctl <- subset(tig.combined,subset=condition=="CTL")

p1<-DimPlot(tig.combined,reduction="pca",group.by = "condition",pt.size=0.01)+
scale_color_manual(values = c("#0072B2","#D55E00"))
p2<-DimPlot(tig.combined.nep,reduction="pca")+
  scale_color_manual(values = c("#0072B2","#D55E00"))
p3<-FeaturePlot(tig.combined.nep,features = "dtd_FLD500",reductio="pca")
  p1+p2+p3
#
# find tig aging marker genes
tig.marker <- FindMarkers(tig.combined.nep,ident.1 = "TIG1-50",ident.2 = "TIG1-20", test.use = "wilcox", min.pct = 0.0,
                          print.bar = TRUE, only.pos = FALSE)
tig.marker$gene <- rownames(tig.marker)
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}
gene.list<-symbol2entrez(rownames(tig.marker[tig.marker$avg_log2FC>0.2 & tig.marker$p_val_adj<0.1,]))
ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))
aging_ego <- ego_result@result[ego_result@result$Description %in% c("aging","cell aging"),]
aging_ego <- ego_result@result[ego_result@result$Description %in% c("epithelial cell proliferation"),]
aging_gene<-unique(unlist(strsplit(aging_ego$geneID,"/")))


