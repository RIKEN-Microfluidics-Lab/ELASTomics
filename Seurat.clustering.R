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
  FeaturePlot(tig,features = c("ITGA8","DLGAP5"), reduction = "umap")
# DLGAP5: tig-1-20
# ITGA8: tig-1-50

tig <- RenameIdents(tig, `0` = "TIG1-50", `1` = "TIG1-20", `2` = "Unknown")
tig <- subset(tig,idents = "Unknown",invert=TRUE)

DimPlot(tig, reduction = "umap",group.by = "condition")+
  DimPlot(tig, reduction = "umap")+
  FeaturePlot(tig,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"), reduction = "umap")
