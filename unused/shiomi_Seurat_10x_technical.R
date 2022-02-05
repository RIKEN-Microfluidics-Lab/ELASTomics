
# show number of counts
VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"),ncol = 2, group.by = "Amu")
tenx <- FetchData(pbmc,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(pbmc,vars="nFeature_RNA")
median(tenx$nFeature_RNA)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 40)
# plot variable features with labels
p1 <- VariableFeaturePlot(pbmc)
p1 <- LabelPoints(plot = p1, points = top10)
p1
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "type")
p2

# running PCA npcs in the RunPCA function must be less than the number of samples
# default is 50

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, npcs=50, features = VariableFeatures(object = pbmc))

#print(pbmc[["pca"]], dims = 1:10, nfeatures = 50)

# gene expression scatter
pca_topcells <- TopCells(object = pbmc[['pca']], balanced = FALSE)
CellScatter(object = pbmc, cell1 = pca_topcells[1], cell2 = pca_topcells[2])


pbmc <- JackStraw(pbmc, num.replicate = 100,dim=40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.3)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
p1 <- DimPlot(pbmc, reduction = "pca", group.by = "ACTB")
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "P53")
p1+p2

#find marker genes in each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.5)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 30)

FeaturePlot(pbmc,reduction = "umap",features="COL1A1", cols = c("blue", "yellow", "red"), max.cutoff = 60, min.cutoff = 0)


rm(p1,p2,pca_topcells,top10,all.genes,tenx)


