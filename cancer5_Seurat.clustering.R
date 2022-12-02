cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer,subset=percent.mt<10)
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

DimPlot(cancer,reduction = "tsne")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CDKN1A","CYBA","percent.mt", "S100A2"),reduction = "tsne", max.cutoff = 4)
cancer <- RenameIdents(cancer, `2` = "unknown", `1` = "MCF10A", `0` = "MCF7")

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

cancer5_merge<- merge(cancer_mcf10,y=cancer_mcf7)
cancer_merge12345 <-merge(cancer_merge1234,y=cancer5_merge)
cancer_merge12345 <- FindVariableFeatures(cancer_merge12345, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge12345)
cancer_merge12345 <- ScaleData(cancer_merge12345, features = all.genes)
cancer_merge12345 <- RunPCA(cancer_merge12345, npcs=20, features = VariableFeatures(object = cancer_merge1234))
DimPlot(cancer_merge12345,group.by = "celltype")


rm(cancer,cancer_mcf10,cancer_mcf7)
rm(cancer5,cancer5_merge)
