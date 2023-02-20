cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer, subset= percent.mt<10 & nCount_RNA>=2000)


#
# PCA and visualize
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=50, features = VariableFeatures(object = cancer))

DimPlot(cancer, reduction = "pca",group.by = "run")+
  DimPlot(cancer,reduction="pca",group.by = "celltype")+
  FeaturePlot(cancer,features = "dtd_FLD004")

cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)

cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.05)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, dims = 1:10)
DimPlot(cancer, reduction = "tsne",group.by = "celltype")+
  DimPlot(subset(cancer,subset=run=="third"),reduction="pca")+
  FeaturePlot(cancer,features = c("dtd_FLD500"), reduction = "umap")+
  FeaturePlot(cancer,features = "percent.mt",reduction = "pca",max.cutoff = 10)

#
# subset third run
cancer3<-subset(cancer,subset=run=="third")
DimPlot(subset(cancer3,subset=percent.mt<5),reduction = "pca")+
  FeaturePlot(subset(cancer3,subset=percent.mt<5),
              features = c("CCDC190","CAVIN3","CD163L1","dtd_FLD004"),reduction = "pca",
              max.cutoff = 4)

cancer3 <- FindVariableFeatures(cancer3, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer3)
cancer3 <- ScaleData(cancer3, features = all.genes)
cancer3 <- RunPCA(cancer3, npcs=50, features = VariableFeatures(object = cancer3))

cancer3 <- FindNeighbors(cancer3, dims = 1:10)
cancer3 <- FindClusters(cancer3, resolution = 0.05)
DimPlot(cancer3, reduction = "umap",group.by = "celltype")+
  FeaturePlot(cancer3,features = "dtd_FLD004")
DimPlot(subset(cancer3,subset=percent.mt<5),reduction = "pca")+
  FeaturePlot(subset(cancer3,subset=percent.mt<5),
              features = c("DUSP1","CYBA","CAVIN3","S100A2"),reduction = "pca",
              max.cutoff = 4)



cancer3 <- RenameIdents(cancer3, `0` = "PC3", `1` = "MCF10A", `2` = "MDAMB231")
#cancer1 <- RenameIdents(cancer1, `0` = "MCF7_NEP", `1` = "PC3_NEP")
cancer_mcf10<-subset(cancer3,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"0V"
cancer_mcf10[["hashtag"]]<-"na"
cancer_mcf10[["condition"]]<-"normal"

cancer_pc3<-subset(cancer3,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"normal"
cancer_pc3[["NEP"]]<-"75V"
cancer_pc3[["hashtag"]]<-"P3NE33"

cancer_mda<-subset(cancer3,idents="MDAMB231")
cancer_mda[["celltype"]]<-"MDAMB231"
cancer_mda[["condition"]]<-"normal"
cancer_mda[["NEP"]]<-"40V"
cancer_mda[["hashtag"]]<-"MM2E32"

cancer3_merge<- merge(cancer_mcf10,y=cancer_pc3)
cancer3_merge<-merge(cancer3_merge,y=cancer_mda)
cancer_merge123 <-merge(cancer_merge,y=cancer3_merge)
cancer_merge123 <- FindVariableFeatures(cancer_merge123, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge123)
cancer_merge123 <- ScaleData(cancer_merge123, features = all.genes)
cancer_merge123 <- RunPCA(cancer_merge123, npcs=20, features = VariableFeatures(object = cancer_merge123))
DimPlot(cancer_merge123,group.by = "celltype")
cancer_merge123 <- RunUMAP(cancer_merge123, dims = 1:10)
cancer_merge123<-RunTSNE(cancer_merge123,dims=1:10)
DimPlot(cancer_merge123,reduction = "tsne",group.by = "celltype")+
  DimPlot(cancer_merge123,reduction = "tsne",group.by = "NEP")+
  FeaturePlot(cancer_merge123,
              features = c("percent.mt"),
              reduction = "tsne",
              max.cutoff = 10)
FeaturePlot(cancer_merge123,
              features = c("CYBA","S100A2","CCDC190","CAVIN3","CD163L1","PLAU","FOSL1","percent.mt"),
              reduction = "tsne",
              max.cutoff = 10)
