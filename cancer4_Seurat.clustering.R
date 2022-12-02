cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer, subset= percent.mt<10)
#
# PCA and visualize
#
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
cancer <- FindClusters(cancer, resolution = 0.1)
#
# Retreiving the results of the preprocessing from the Seurat object
#
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(cancer, reduction = "pca",group.by = "celltype")+
  DimPlot(subset(cancer,subset=run=="forth"),reduction="pca")+
  FeaturePlot(cancer,features = c("dtd_FLD500"), reduction = "umap")
#
# subset forth run
#
cancer4<-subset(cancer,subset=run=="forth")
cancer4 <- FindVariableFeatures(cancer4, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer4)
cancer4 <- ScaleData(cancer4, features = all.genes)
cancer4 <- RunPCA(cancer4, npcs=50, features = VariableFeatures(object = cancer4))
cancer4 <- FindNeighbors(cancer4, dims = 1:10)
cancer4 <- FindClusters(cancer4, resolution = 0.04)
DimPlot(cancer4, reduction = "pca")+
  FeaturePlot(cancer4,features = "dtd_FLD004")
DimPlot(cancer4,reduction = "pca")+
  FeaturePlot(cancer4, features = c("CAVIN3","S100A2","DUSP1"),reduction = "pca", max.cutoff = 4)
cancer4 <- RenameIdents(cancer4, `1` = "PC3", `3` = "PC3", `0` = "MCF10A", `2` = "MDAMB231")

cancer_mcf10<-subset(cancer4,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"40V"
cancer_mcf10[["hashtag"]]<-"M1AE31"
cancer_mcf10[["condition"]]<-"normal"

cancer_pc3<-subset(cancer4,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"mbcd"
cancer_pc3[["NEP"]]<-"40V"
cancer_pc3[["hashtag"]]<-"P3ME34"

cancer_mda<-subset(cancer4,idents="MDAMB231")
cancer_mda[["celltype"]]<-"MDAMB231"
cancer_mda[["condition"]]<-"normal"
cancer_mda[["NEP"]]<-"0V"
cancer_mda[["hashtag"]]<-"na"

cancer4_merge<- merge(cancer_mcf10,y=cancer_pc3)
cancer4_merge<-merge(cancer4_merge,y=cancer_mda)
cancer_merge1234 <-merge(cancer_merge123,y=cancer4_merge)
cancer_merge1234 <- FindVariableFeatures(cancer_merge1234, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge1234)
cancer_merge1234 <- ScaleData(cancer_merge1234, features = all.genes)
cancer_merge1234 <- RunPCA(cancer_merge1234, npcs=20, features = VariableFeatures(object = cancer_merge1234))
DimPlot(cancer_merge1234,group.by = "celltype")

rm(cancer,cancer_mcf10,cancer_mcf7,cancer_mda,cancer_merge,cancer_merge123,cancer_pc3)
rm(cancer1,cancer2,cancer2_mcf7,cancer2_pc3,cancer3_merge,cancer3,cancer4,cancer4_merge)