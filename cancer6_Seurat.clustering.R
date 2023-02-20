cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer,subset=percent.mt<10 & nCount_RNA>=2000)
#
# PCA and visualize
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=50, features = VariableFeatures(object = cancer))

DimPlot(cancer, reduction = "pca",group.by = "run")+
  DimPlot(cancer,reduction="pca",group.by = "celltype")+
  DimPlot(subset(cancer,subset=run==c("fifth","sixth")),reduction="pca")+
  FeaturePlot(cancer,features = "dtd_FLD004")

cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)

cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.1)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(cancer, reduction = "pca",group.by = "celltype")+
  DimPlot(subset(cancer,subset=run=="forth"),reduction="pca")+
  FeaturePlot(cancer,features = c("dtd_FLD500"), reduction = "umap")

#
# subset third run
cancer6<-subset(cancer,subset=run=="sixth")

cancer6 <- FindVariableFeatures(cancer6, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer6)
cancer6 <- ScaleData(cancer6, features = all.genes)
cancer6 <- RunPCA(cancer6, npcs=50, features = VariableFeatures(object = cancer6))

cancer6 <- FindNeighbors(cancer6, dims = 1:10)
cancer6 <- FindClusters(cancer6, resolution = 0.04)
DimPlot(cancer6, reduction = "pca")+
  FeaturePlot(cancer6,features = "dtd_FLD004")
DimPlot(cancer6,reduction = "pca")+
  FeaturePlot(cancer6,
              features = c("DUSP1","CYBA","dtd_FLD004"),reduction = "pca",
              max.cutoff = 4)
DimPlot(cancer,reduction = "pca",group.by = "celltype")+
  FeaturePlot(cancer,
              features = c("DUSP1","CYBA"),reduction = "pca",
              max.cutoff = 4)


cancer6 <- RenameIdents(cancer6,`0` = "MCF10A", `1` = "MCF7")

#cancer1 <- RenameIdents(cancer1, `0` = "MCF7_NEP", `1` = "PC3_NEP")
cancer_mcf10<-subset(cancer6,idents = "MCF10A")
cancer_mcf10[["celltype"]]<-"MCF10A"
cancer_mcf10[["NEP"]]<-"40V"
#cancer_mcf10[["hashtag"]]<-"M1AE31"
cancer_mcf10[["condition"]]<-"cytD"

#cancer_pc3<-subset(cancer4,idents="PC3")
#cancer_pc3[["celltype"]]<-"PC3"
#cancer_pc3[["condition"]]<-"mbcd"
#cancer_pc3[["NEP"]]<-"40V"
#cancer_pc3[["hashtag"]]<-"P3ME34"

cancer_mcf7<-subset(cancer6,idents="MCF7")
cancer_mcf7[["celltype"]]<-"MCF7"
cancer_mcf7[["condition"]]<-"normal"
cancer_mcf7[["NEP"]]<-"40V"
cancer_mcf7[["hashtag"]]<-"M7AE37"

cancer6_merge<- merge(cancer_mcf10,y=cancer_mcf7)
#cancer4_merge<-merge(cancer4_merge,y=cancer_mda)
cancer_merge123456 <-merge(cancer_merge12345,y=cancer6_merge)
cancer_merge123456 <- FindVariableFeatures(cancer_merge123456, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge123456)
cancer_merge123456 <- ScaleData(cancer_merge123456, features = all.genes)
cancer_merge123456 <- RunPCA(cancer_merge123456, npcs=20, features = VariableFeatures(object = cancer_merge123456))
DimPlot(cancer_merge123456,group.by = "celltype")
cancer_merge123456 <- RunUMAP(cancer_merge123456, dims = 1:10)
cancer_merge123456 <- RunTSNE(cancer_merge123456,dims=1:10)

  FeaturePlot(subset(subset(cancer_merge123456,subset=condition=="normal"),subset=NEP=="75V",invert=TRUE),
              features = c("CYBA","S100A2","CCDC190","CAVIN3","CD163L1","PLAU","FOSL1","percent.mt"),
              reduction = "tsne",
              max.cutoff = 10)+
    DimPlot(subset(subset(cancer_merge123456,subset=condition=="normal"),subset=NEP=="75V",invert=TRUE),
            reduction = "tsne",group.by = "celltype")+
    DimPlot(subset(subset(cancer_merge123456,subset=condition=="normal"),subset=NEP=="75V",invert=TRUE),
            reduction = "tsne",group.by = "NEP")
  
  
  FeaturePlot(subset(cancer_merge123456,subset=condition=="normal"),
              features = c("CYBA","S100A2","CCDC190","CAVIN3","PLAU","FOSL1"),
              reduction = "tsne",
              max.cutoff = 10)

  p1<-FeaturePlot(subset(subset(cancer_merge123456,subset=NEP== "0V"),subset=condition=="normal"),
                features = c("percent.mt"),
                reduction = "tsne",
                max.cutoff = 10)
    p2<-FeaturePlot(subset(subset(cancer_merge123456,subset=NEP== "0V",invert=TRUE),subset=condition=="normal"),
                features = c("percent.mt"),
                reduction = "tsne",
                max.cutoff = 10)
    
    p3<-DimPlot(subset(cancer_merge123456,subset=percent.mt<10& nCount_RNA>=2000 & condition=="normal"),reduction = "tsne",group.by = "celltype")
    p4<-DimPlot(subset(cancer_merge123456,subset=percent.mt<5& nCount_RNA>=2000 & condition=="normal"),reduction = "tsne",group.by = "NEP")
    p5<-DimPlot(subset(cancer_merge123456,subset=percent.mt<10& nCount_RNA>=2000 & condition=="normal"),reduction = "tsne",group.by = "NEP")
    p6<-DimPlot(subset(cancer_merge123456,subset=percent.mt<5& nCount_RNA>=2000 & condition=="normal"),reduction = "tsne",group.by = "NEP")
  gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3)
  
      FeaturePlot(subset(subset(cancer_merge123456,subset=percent.mt<5),subset=condition=="normal"),features ="nCount_RNA",
              reduction = "tsne",max.cutoff = 10000)
  
DimPlot(subset(cancer_merge123456,subset=percent.mt<5),reduction = "umap",group.by = "celltype")+
  FeaturePlot(subset(cancer_merge12345,subset=percent.mt<5),
              features = c("CYBA","S100A2","CCDC190","CAVIN3","CD163L1","PLAU","FOSL1"),
              reduction = "umap",
              max.cutoff = 4)

rm(cancer,cancer_mcf10,cancer_mcf7,cancer_merge12346)
rm(cancer6,cancer6_merge)
rm(p1,p2,p3,p4,p5,p6)
