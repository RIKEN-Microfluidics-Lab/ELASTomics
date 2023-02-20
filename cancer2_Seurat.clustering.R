cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer, subset= percent.mt<10)
#
# merge all the tig data
#
cancer <- NormalizeData(cancer, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
#
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)
#
# subset second run
#
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.1)
cluster = as.numeric(Idents(cancer))
FeatureScatter(cancer,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA")
DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CDKN1A","CYBA","S100A2","dtd_FLD004","dtd_P3NE35","dtd_MC7E37", "nFeature_RNA","percent.mt"),reduction = "umap")
cancer <- RenameIdents(cancer, `0` = "dead", `1` = "PC3", `2` = "MCF7", `3` = "PC3")
cancer[["celltype"]]<-Idents(cancer)

cancer2_mcf7<-subset(cancer,idents = "MCF7")
cancer2_mcf7[["celltype"]]<-"MCF7"
cancer2_mcf7[["NEP"]]<-"0V"
cancer2_mcf7[["hashtag"]]<-"na"
cancer2_mcf7[["condition"]]<-"normal"

cancer2_pc3<-subset(cancer,idents="PC3")
cancer2_pc3[["celltype"]]<-"PC3"
cancer2_pc3[["NEP"]]<-"40V"
cancer2_pc3[["hashtag"]]<-"36"
cancer2_pc3[["condition"]]<-"CytoD"
cancer2<-merge(cancer2_mcf7,y=cancer2_pc3)
rm(cancer2_mcf7, cancer2_pc3)