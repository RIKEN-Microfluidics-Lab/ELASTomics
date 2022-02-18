#https://satijalab.org/seurat/articles/integration_introduction.html
tig.list <- SplitObject(tig, split.by = "condition")
tig.list <- lapply(X = tig.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 400)
})
features <- SelectIntegrationFeatures(object.list = tig.list)

tig.anchors <- FindIntegrationAnchors(object.list = tig.list, anchor.features = features)
tig.combined <- IntegrateData(anchorset = tig.anchors)
DefaultAssay(tig.combined) <- "integrated"
tig.combined <- ScaleData(tig.combined, verbose = FALSE)
tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindClusters(tig.combined, resolution = 0.5)
p1 <- DimPlot(tig.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(tig.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
# remove low quality samples
tig.combined <-subset(tig.combined,idents=c(6),invert=TRUE)
tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindClusters(tig.combined, resolution = 0.1)
p1 <- DimPlot(tig.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(tig.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

tig.combined <- NormalizeData(tig.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
FeaturePlot(tig.combined,features=c("dtd_FLD004","CDKN1A"))
tig.combined.nep <- subset(tig.combined,subset=condition=="NEP")


tig.combined <- RenameIdents(tig.combined, `0` = "TIG1-50", `1` = "TIG1-20")
FeaturePlot(tig.combined,features=c("dtd_FLD004","CDKN1A"))
DimPlot(tig.combined)
tig.50 <- subset(tig.combined,idents="TIG1-50")
tig.50<- CellCycleScoring(tig.50,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)

tig.20 <- subset(tig.combined,idents="TIG1-20")
tig.20<- CellCycleScoring(tig.20,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
tig.50 <- NormalizeData(tig.50,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig.20 <- NormalizeData(tig.20,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
DimPlot(tig.20,reduction="umap")+FeaturePlot(tig.20,features = "dtd_FLD004")+
  DimPlot(tig.50,reduction="umap")+FeaturePlot(tig.50,features = "dtd_FLD004")

