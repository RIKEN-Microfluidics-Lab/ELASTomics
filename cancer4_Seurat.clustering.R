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
# subset forth run
#
cancer <- RunUMAP(cancer, dims = 1:10)
cancer <- RunTSNE(cancer, npcs=50, features = VariableFeatures(object = cancer))
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.05)
cluster = as.numeric(Idents(cancer))
DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CYBA","S100A2","dtd_FLD004","dtd_MM2E32", "dtd_P3HE33", "nFeature_RNA","percent.mt"),reduction = "umap")
cancer <- RenameIdents(cancer, `0` = "MCF10A", `1` = "PC3", `2` = "MDAMB231", `3` = "PC3")

cancer_mcf10<-subset(cancer,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"40V"
cancer_mcf10[["hashtag"]]<-"M1AE31"
cancer_mcf10[["condition"]]<-"normal"

cancer_pc3<-subset(cancer,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"mbcd"
cancer_pc3[["NEP"]]<-"40V"
cancer_pc3[["hashtag"]]<-"P3ME34"

cancer_mda<-subset(cancer,idents="MDAMB231")
cancer_mda[["celltype"]]<-"MDAMB231"
cancer_mda[["condition"]]<-"normal"
cancer_mda[["NEP"]]<-"0V"
cancer_mda[["hashtag"]]<-"na"

cancer4 <- merge(cancer_mcf10,y=cancer_pc3)
cancer4 <- merge(cancer4,y=cancer_mda)
rm(cancer, cancer_pc3, cancer_mda, cancer_mcf10)