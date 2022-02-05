
tig <- JackStraw(tig, num.replicate = 100)
tig <- ScoreJackStraw(tig, dims = 1:20)
JackStrawPlot(tig, dims = 1:20)
ElbowPlot(tig)

tig <- FindNeighbors(tig, dims = 1:19)
tig <- FindClusters(tig, resolution = 0.1)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(tig))
tig <- RunUMAP(tig, dims = 1:19)
DimPlot(tig, reduction = "umap",group.by = "condition")+
  DimPlot(tig, reduction = "umap")+
  FeaturePlot(tig,features = "DLGAP5", reduction = "umap")
#  FeaturePlot(tig,features = "percent.mt", reduction = "umap")


#find marker genes in each cluster
tig.markers <- FindAllMarkers(tig, only.pos = FALSE )
