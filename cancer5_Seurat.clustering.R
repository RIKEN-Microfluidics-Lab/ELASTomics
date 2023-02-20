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
DimPlot(cancer, reduction = "pca",group.by = "run")
cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)
#
# subset fifth run
#
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.1)
cluster = as.numeric(Idents(cancer))

DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CYBA","S100A2","dtd_FLD004", "nFeature_RNA","percent.mt"),reduction = "umap")
DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("dtd_FLD004"),reduction = "umap")

cancer <- RenameIdents(cancer, `0` = "MCF10A", `1` = "MCF7", `2` = "MCF7", `3` = "dead")

cancer_mcf10<-subset(cancer,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"40V"
cancer_mcf10[["hashtag"]]<-"MCF10A"
cancer_mcf10[["condition"]]<-"normal"

cancer_mcf7<-subset(cancer,idents="MCF7")
cancer_mcf7[["celltype"]]<-"MCF7"
cancer_mcf7[["condition"]]<-"normal"
cancer_mcf7[["NEP"]]<-"0V"
cancer_mcf7[["hashtag"]]<-"na"

cancer5<- merge(cancer_mcf10,y=cancer_mcf7)
rm(cancer,cancer_mcf10,cancer_mcf7)

