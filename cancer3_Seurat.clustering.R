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
# subset third run
#
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.1)
cluster = as.numeric(Idents(cancer))
FeatureScatter(cancer,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA")
DimPlot(cancer,reduction = "tsne")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CYBA","S100A2","dtd_FLD004","dtd_MM2E32", "dtd_P3HE33", "nFeature_RNA","percent.mt"),reduction = "tsne")
cancer <- RenameIdents(cancer, `0` = "PC3", `1` = "MCF10A", `2` = "MDAMB231", `3` = "PC3", `4` = "MCF10A")

cancer_mcf10<-subset(cancer,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"0V"
cancer_mcf10[["hashtag"]]<-"na"
cancer_mcf10[["condition"]]<-"normal"

cancer_pc3<-subset(cancer,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"normal"
cancer_pc3[["NEP"]]<-"75V"
cancer_pc3[["hashtag"]]<-"P3NE33"

cancer_mda<-subset(cancer,idents="MDAMB231")
cancer_mda[["celltype"]]<-"MDAMB231"
cancer_mda[["condition"]]<-"normal"
cancer_mda[["NEP"]]<-"40V"
cancer_mda[["hashtag"]]<-"MM2E32"

cancer3 <- merge(cancer_mcf10,y=cancer_pc3)
cancer3 <- merge(cancer3,y=cancer_mda)
rm(cancer, cancer_pc3, cancer_mda, cancer_mcf10)