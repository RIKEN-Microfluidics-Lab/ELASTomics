cancer_merge<-subset(subset(subset(subset(cancer_merge1234,subset=percent.mt<10),subset=nCount_RNA>=2000),
                     subset=NEP=="75V",invert=TRUE),subset=condition=="normal")

cancer_merge <- FindVariableFeatures(cancer_merge, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge)
cancer_merge <- ScaleData(cancer_merge, features = all.genes)
cancer_merge <- RunPCA(cancer_merge, npcs=20, features = VariableFeatures(object = cancer_merge))
DimPlot(cancer_merge,group.by = "celltype")
cancer_merge <- RunUMAP(cancer_merge, dims = 1:10)
cancer_merge <- RunTSNE(cancer_merge,dims=1:10)

FeaturePlot(cancer_merge,
            features = c("CYBA","S100A2","CCDC190","CAVIN3","CD163L1","PLAU","FOSL1","percent.mt"),
            reduction = "tsne",
            max.cutoff = 10)+
  DimPlot(cancer_merge,reduction = "tsne",group.by = "celltype")+
  DimPlot(cancer_merge,reduction = "tsne",group.by = "NEP")


FeaturePlot(subset(cancer_merge,subset=percent.mt<10),
            features = c("CYBA","S100A2","CCDC190","CAVIN3","PLAU","FOSL1"),
            reduction = "tsne",
            max.cutoff = 10)

p1<-FeaturePlot(subset(subset(cancer_merge,subset=NEP== "0V"),subset=condition=="normal"),
                features = c("percent.mt"),
                reduction = "tsne",
                max.cutoff = 10)
p2<-FeaturePlot(subset(subset(cancer_merge,subset=NEP== "0V",invert=TRUE),subset=condition=="normal"),
                features = c("percent.mt"),
                reduction = "tsne",
                max.cutoff = 10)
p3<-DimPlot(subset(subset(cancer_merge,subset=percent.mt<10),subset=condition=="normal"),reduction = "tsne",group.by = "celltype")
p4<-DimPlot(subset(cancer_merge,subset=condition=="normal"),reduction = "tsne",group.by = "NEP")
p5<-DimPlot(subset(subset(cancer_merge,subset=percent.mt<10),subset=condition=="normal"),reduction = "tsne",group.by = "NEP")
p6<-DimPlot(subset(subset(cancer_merge,subset=percent.mt<5),subset=condition=="normal"),reduction = "tsne",group.by = "NEP")
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3)

FeaturePlot(subset(cancer_merge,subset=percent.mt<5 & condition=="normal"),features ="nCount_RNA",
            reduction = "tsne",max.cutoff = 10000)


  FeaturePlot(subset(cancer_merge,subset=percent.mt<5 & condition=="normal" & NEP=="40V"),
              features = c("CYBA","S100A2","CCDC190","CAVIN3","CD163L1","PLAU","FOSL1","percent.mt"),
              reduction = "tsne",
              max.cutoff = 4)+
    DimPlot(subset(cancer_merge,subset=percent.mt<5 & condition=="normal" & NEP=="40V"),
            reduction = "tsne",group.by = "celltype")+
    DimPlot(subset(cancer_merge,subset=percent.mt<5 & condition=="normal"),
            reduction = "tsne",group.by = "NEP")

