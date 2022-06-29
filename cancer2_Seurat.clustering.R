cancer <- NormalizeData(cancer, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
#cancer <- subset(cancer, subset= percent.mt<5)
#
#
# merge all the tig data
#
cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(cancer), 40)
# plot variable features with labels
#plot1 <- VariableFeaturePlot(cancer)
#plot1 <- LabelPoints(plot = cancer, points = top10)
#plot1
# normalize the dtd data with "RC" option.
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e6)
#
# PCA and visualize
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=20, features = VariableFeatures(object = cancer))
DimPlot(cancer, reduction = "pca",group.by = "run")

cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)

cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.02)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(cancer,reduction="pca")+
  FeaturePlot(cancer,features = c("dtd_P3CE36","dtd_MC7E37","dtd_FLD004"),
              reduction = "pca",slot="counts",
              max.cutoff = 100)
# DLGAP5: cancer-1-20
# ITGA8: cancer-1-50

cancer <- RenameIdents(cancer, `1` = "PC3_NEP", `0` = "MCF7_CTL")
#FindAllMarkers(cancer)
cancer2_mcf7<-subset(cancer,idents = "MCF7_CTL")
cancer2_mcf7[["celltype"]]<-"MCF7"
cancer2_mcf7[["NEP"]]<-"0V"
cancer2_mcf7[["hashtag"]]<-"na"
cancer2_mcf7[["condition"]]<-"normal"
cancer2_pc3<-subset(cancer,idents="PC3_NEP")
cancer2_pc3[["celltype"]]<-"PC3"
cancer2_pc3[["NEP"]]<-"40V"
cancer2_pc3[["hashtag"]]<-"36"
cancer2_pc3[["condition"]]<-"CytoD"
cancer2<-merge(cancer2_mcf7,y=cancer2_pc3)
#cancer2 <- subset(cancer2, subset=percent.mt <10)
cancer2 <- FindVariableFeatures(cancer2, selection.method = "vst", nfeatures = 100)
#
# PCA and visualize
all.genes <- rownames(cancer2)
cancer2 <- ScaleData(cancer2, features = all.genes)
cancer2 <- RunPCA(cancer2, npcs=20, features = VariableFeatures(object = cancer2))
DimPlot(cancer2, reduction = "pca",group.by = "celltype")

#DimPlot(cancer2, reduction = "umap")+
  DimPlot(subset(cancer2, subset= percent.mt<5),reduction="pca")+
    FeaturePlot(subset(cancer2, subset= percent.mt<5),features =c("CYBA","SERPINA3","PLAU","FOSL1"), reduction = "pca")+
    VlnPlot(subset(cancer2, subset= percent.mt<5),features = "dtd_FLD004", group.by = "celltype")+
    VlnPlot(subset(cancer2, subset= percent.mt<5),features = "FOSL1", group.by = "celltype")
FeatureScatter(cancer2,feature1 = "dtd_P3CE36",feature2 = "dtd_FLD004",slot = "counts")
VlnPlot(cancer2,features =c("percent.mt","nCount_RNA"), group.by = "celltype")
