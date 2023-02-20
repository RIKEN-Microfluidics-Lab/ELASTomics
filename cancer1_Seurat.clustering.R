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
# subset first run
#
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.1)
cluster = as.numeric(Idents(cancer))
DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("DUSP1","CDKN1A","CYBA","S100A2","dtd_FLD004","dtd_P3NE35","dtd_MC7E37", "nFeature_RNA","percent.mt"),reduction = "umap")
cancer <- RenameIdents(cancer,`0` = "PC3", `1` = "PC3", `2` = "MCF7", `3` = "merge")
cancer[["celltype"]]<-Idents(cancer)

cancer_mcf7<-subset(cancer,idents = "MCF7")
cancer_mcf7[["celltype"]]<-"MCF7"
cancer_mcf7[["NEP"]]<-"40V"
cancer_mcf7[["hashtag"]]<-"MC7E37"
cancer_mcf7[["condition"]]<-"normal"
threshold <-colMeans(FetchData(cancer_mcf7,"dtd_P3NE35"))

cancer_pc3<-subset(cancer,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"normal"
cancer_pc3[["NEP"]]<-"0V"
nep<-data.frame(cancer_pc3[["NEP"]])
nep[colnames(subset(cancer_pc3,subset=dtd_P3NE35>threshold)),] <-"40V"
cancer_pc3[["NEP"]]<-nep$NEP
cancer_pc3[["hashtag"]]<-"P3NE35"


cancer1<-merge(cancer_mcf7,y=cancer_pc3)
rm(cancer_mcf7, cancer_pc3, cancer)
