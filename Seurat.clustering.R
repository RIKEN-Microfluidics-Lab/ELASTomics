cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=20, features = VariableFeatures(object = cancer))
DimPlot(cancer, reduction = "pca",group.by = "run")+
  DimPlot(cancer,reduction="pca",group.by = "celltype")+
  DimPlot(subset(cancer,subset=NEP=="0V"),reduction="pca",group.by = "celltype")+
  FeaturePlot(subset(cancer,subset=condition=="normal"),features = "dtd_FLD004")

cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)

cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.02)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(cancer, reduction = "pca")+
  DimPlot(subset(cancer,subset=run=="second"),reduction="pca",group.by = "celltype")+
  FeaturePlot(cancer,features = c("dtd_FLD500"), reduction = "pca")
# DLGAP5: cancer-1-20
# ITGA8: cancer-1-50

cancer <- RenameIdents(cancer, `0` = "MCF7", `1` = "PC3")
dtd.mtx <-data.frame(t(cancer[["DTD"]]@counts))
dtd.mtx$celltype <-"MCF7"
dtd.mtx[dtd.mtx$P3NE35>dtd.mtx$MC7E37 | dtd.mtx$MC7E37==0,]$celltype <-"PC3"
cancer[["tag"]]<-dtd.mtx$celltype
DimPlot(subset(cancer,subset=run=="first"), group.by = "celltype",reduction="pca")+
FeatureScatter(subset(cancer,subset=run=="first"),feature1 = "dtd_MC7E37",feature2="dtd_P3NE35",group.by = "celltype")+
  FeatureScatter(subset(cancer,subset=run=="first"),feature1 = "dtd_MC7E37",feature2="dtd_P3NE35",slot="counts",group.by = "celltype")
#cancer <- subset(cancer,idents = "Unknown",invert=TRUE)
FeatureScatter(subset(cancer,subset=run=="first"),
            feature1 = "dtd_MC7E37", feature2 = "dtd_P3NE35",
            pt.size = 0.1,slot = "counts")+scale_x_log10()+scale_y_log10()
FeaturePlot(subset(cancer,subset=run=="first"),
            features = "dtd_MC7E37", reduction = "pca",
            max.cutoff = 4,min.cutoff = 2)+
  FeaturePlot(subset(cancer,subset=run=="first"),
              features = "dtd_P3NE35", reduction = "pca",
              max.cutoff = 4,min.cutoff = 0)  
VlnPlot(subset(cancer,subset=run=="first"),features = "dtd_FLD500",
        group.by = "celltype")

DimPlot(cancer, reduction = "umap",group.by = "condition")+
DimPlot(cancer, reduction = "umap")+
FeaturePlot(cancer,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"), reduction = "umap")

#  FeaturePlot(cancer,features = "percent.mt", reduction = "umap")
#find marker genes in each cluster
#cancer.markers <- FindAllMarkers(cancer, only.pos = FALSE )
