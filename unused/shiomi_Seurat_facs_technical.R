
# show number of counts
VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"),
        ncol = 2,group.by = "cell")
VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"),
        ncol = 2,group.by = "gate")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "gate" )
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "cell" )


count_summary <- pbmc[[c("nCount_RNA","nFeature_RNA","gate","cell")]]
count_summary_mean <- count_summary %>%
  dplyr::group_by(gate,cell) %>%
  dplyr::summarise(nCount.mean=mean(nCount_RNA), nFeature.mean=mean(nFeature_RNA))
# number of UMI counts under various conditions
ggplot(count_summary_mean,aes(x=factor(gate),y=factor(cell),
                              color=nCount.mean,size=nCount.mean))+
  geom_point()+
  scale_size(range = c(10, 20))

tenx <- FetchData(pbmc,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(pbmc,vars="nFeature_RNA")
median(tenx$nFeature_RNA)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 40)
# plot variable features with labels
plot1 <- VariableFeaturePlot(pbmc)
plot1 <- LabelPoints(plot = plot1, points = top10)
plot1
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" )
plot2

# running PCA npcs in the RunPCA function must be less than the number of samples
# default is 50

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, npcs=10, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:2, nfeatures = 50)

# plot plates/dishes/gates/pools/rtid
p1 <- DimPlot(pbmc, reduction = "pca",group.by = "plates")
p2 <- DimPlot(pbmc, reduction = "pca",group.by = "cell")
p3 <- DimPlot(pbmc, reduction = "pca",group.by = "gate")
p1+p2+p3

DimPlot(pbmc, reduction = "pca",group.by = "pool")
DimPlot(pbmc, reduction = "pca",group.by = "rtid")

# gene expression scatter
pca_topcells <- TopCells(object = pbmc[['pca']], balanced = FALSE)
CellScatter(object = pbmc, cell1 = pca_topcells[1], cell2 = pca_topcells[2])


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:10)
JackStrawPlot(pbmc, dims = 1:10)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.3)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:7)
p1 <- DimPlot(pbmc, reduction = "pca")
p2 <- DimPlot(pbmc, reduction = "umap")
p1+p2

#find marker genes in each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5)

#pbmc.markers <- FindMarkers(pbmc,logfc.threshold = 0.2,ident.1=colnames(subset(pbmc,subset=gate=="RG")))
#pbmc.markers <- FindMarkers(pbmc,logfc.threshold = 0.2,ident.1=colnames(subset(pbmc,subset=cell=="HEA")))

#features <-"Venus"
#p1 <- DimPlot(pbmc, reduction = "pca",group.by = "dish")
#p2 <- DimPlot(pbmc, reduction = "pca",group.by = "gate")
#p3 <- FeaturePlot(pbmc, features = features,reduction = "pca")
#p1+p2+p3

#p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ,group.by = "gate")
#p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ,group.by = "cell")
#p1+p2
#VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"), ncol = 2,group.by = "gate")
#VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"), ncol = 2,group.by = "cell")

rm(p1,p2,p3,plot1,plot2,count_summary,count_summary_mean,pca_topcells,top10,all.genes,tenx)


