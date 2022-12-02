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
DimPlot(cancer, reduction = "pca",group.by = "run")+
  FeaturePlot(cancer,features = "dtd_FLD004")
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
cancer <- FindClusters(cancer, resolution = 0.04)
cluster = as.numeric(Idents(cancer))

DimPlot(cancer,reduction = "umap")+
  FeaturePlot(cancer, features = c("DUSP1","CAVIN3","CDKN1A","CYBA","S100A2", "nCount_RNA"),reduction = "umap")
cancer <- RenameIdents(cancer, `2` = "MCF10A", `1` = "MCF7", `0` = "unknown")

cancer_mcf10<-subset(cancer,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"40V"
cancer_mcf10[["hashtag"]]<-"MCF10C"
cancer_mcf10[["condition"]]<-"CytoD"

cancer_mcf7<-subset(cancer,idents="MCF7")
cancer_mcf7[["celltype"]]<-"MCF7"
cancer_mcf7[["condition"]]<-"normal"
cancer_mcf7[["NEP"]]<-"40V"
cancer_mcf7[["hashtag"]]<-"MC7E37"

cancer6_merge<- merge(cancer_mcf10,y=cancer_mcf7)

cancer_merge123456 <-merge(cancer_merge12345,y=cancer6_merge)
cancer_merge123456 <- subset(cancer_merge123456,subset=nFeature_RNA>2500)
cancer_merge123456 <- FindVariableFeatures(cancer_merge123456, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge123456)
cancer_merge123456 <- ScaleData(cancer_merge123456, features = all.genes)
cancer_merge123456 <- RunPCA(cancer_merge123456, npcs=20, features = VariableFeatures(object = cancer_merge123456))
DimPlot(cancer_merge123456,group.by = "celltype")
cancer_merge123456 <- RunUMAP(cancer_merge123456, dims = 1:10)
cancer_merge123456 <- RunTSNE(cancer_merge123456,dims=1:10)

  FeaturePlot(cancer_merge123456,
              features = c("CYBA","S100A2","CAVIN3","FOSL1","CDKN1A","DUSP1"),
              reduction = "tsne",
              max.cutoff = 10)+
    DimPlot(cancer_merge123456,reduction = "tsne",group.by = "celltype")+
    DimPlot(cancer_merge123456,reduction = "tsne",group.by = "NEP")

rm(cancer,cancer_mcf10,cancer_mcf7,cancer6_merge, cancer6, cancer_merge1234, cancer_merge12345)
