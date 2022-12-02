cancer <- NormalizeData(cancer, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
cancer <- subset(cancer, subset= percent.mt<10)
#
# merge all the tig data
#
cancer <- FindVariableFeatures(cancer, selection.method = "vst", nfeatures = 300)
cancer <- NormalizeData(cancer,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
#
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, npcs=20, features = VariableFeatures(object = cancer))
DimPlot(subset(cancer,subset=run=="first"), reduction = "pca")+
  DimPlot(subset(cancer,subset=run=="second"),reduction="pca",group.by = "celltype")+
  DimPlot(cancer,reduction="pca")+
  DimPlot(cancer,reduction="pca",group.by = "celltype")
FeaturePlot(subset(cancer,subset=run=="first"), reduction = "pca",features = "percent.mt")+
  FeaturePlot(subset(cancer,subset=run=="second"),reduction="pca",features="percent.mt")
cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:20)
ElbowPlot(cancer)
cancer2<-subset(cancer,subset=run=="second")
cancer <- FindNeighbors(cancer, dims = 1:10)
cancer <- FindClusters(cancer, resolution = 0.03)
#
# Retreiving the results of the preprocessing from the Seurat object
#
cluster = as.numeric(Idents(cancer))
cancer <- RunUMAP(cancer, dims = 1:10)
DimPlot(subset(cancer,subset=run=="first"), reduction = "pca")+
  DimPlot(subset(cancer,subset=run=="second"),reduction="pca",group.by = "celltype")+
  DimPlot(cancer, reduction = "pca")
DimPlot(cancer,reduction="pca")+
  VlnPlot(subset(cancer,subset=run=="first"),features = c("dtd_P3NE35","dtd_MC7E37","CYBA","S100A2","PLAU","FOSL1"))
FeaturePlot(subset(cancer,subset=run=="first"),reduction="pca",
            features=c("dtd_P3NE35","dtd_MC7E37","CYBA","S100A2","PLAU","FOSL1"),
            slot="counts")
cancer1<-subset(cancer,subset=run=="first")
cancer1 <- RenameIdents(cancer1,`0` = "PC3", `1` = "PC3", `2` = "MCF7")
cancer1[["celltype"]]<-Idents(cancer1)

cancer_mcf7<-subset(cancer1,idents = "MCF7")
cancer_mcf7[["celltype"]]<-"MCF7"
cancer_mcf7[["NEP"]]<-"40V"
cancer_mcf7[["hashtag"]]<-"MC7E37"
cancer_mcf7[["condition"]]<-"normal"
threshold <-colMeans(FetchData(cancer_mcf7,"dtd_P3NE35"))

cancer_pc3<-subset(cancer1,idents="PC3")
cancer_pc3[["celltype"]]<-"PC3"
cancer_pc3[["condition"]]<-"normal"
cancer_pc3[["NEP"]]<-"0V"
nep<-data.frame(cancer_pc3[["NEP"]])
nep[colnames(subset(cancer_pc3,subset=dtd_P3NE35>threshold)),] <-"40V"
cancer_pc3[["NEP"]]<-nep$NEP
cancer_pc3[["hashtag"]]<-"P3NE35"

#FeaturePlot(cancer_pc3,features="dtd_P3NE35",reduction="pca",slot="counts")
#cancer_pc3_NEP <- subset(cancer_pc3)
#cancer_pc3_NEP[["NEP"]]<-"40V"
#cancer_pc3_NEP[["hashtag"]]<-"P3NE35"
#cancer_pc3_CTL <- subset(cancer_pc3,subset=dtd_P3NE35==0)
#cancer_pc3_NEP[["NEP"]]<-"0V"
#cancer_pc3_NEP[["hashtag"]]<-"NA"

cancer1<-merge(cancer_mcf7,y=cancer_pc3)
cancer_merge <- merge(cancer1,y=cancer2)
cancer_merge[["celltype_NEP_condition"]]<-paste0(unlist(cancer_merge[["celltype"]]),"_",unlist(cancer_merge[["NEP"]]),"_",unlist(cancer_merge[["condition"]]))
all.genes <- rownames(cancer_merge)
cancer_merge <- ScaleData(cancer_merge, features = all.genes)
cancer_merge <- RunPCA(cancer_merge, npcs=20, features = VariableFeatures(object = cancer))
cancer_merge <- NormalizeData(cancer_merge,assay="DTD",normalization.method = "CLR",scale.factor = 1e6)
DimPlot(subset(cancer_merge,subset=run=="first"), reduction = "pca",group.by = "celltype_NEP_condition")+
  FeaturePlot(subset(cancer_merge,subset=run=="first"),reduction="pca",features=c("dtd_P3NE35","dtd_MC7E37","dtd_P3CE36"))
DimPlot(subset(cancer_merge,subset=run=="second"), reduction = "pca",group.by = "celltype_NEP_condition")+
    FeaturePlot(subset(cancer_merge,subset=run=="second"),reduction="pca",features=c("dtd_P3NE35","dtd_MC7E37","dtd_P3CE36"))

FeaturePlot(cancer_merge,reduction="pca",features=c("dtd_P3NE35","dtd_MC7E37","dtd_P3CE36"))
VlnPlot(cancer_merge,features=c("dtd_FLD004","dtd_FLD500","percent.mt"),group.by = "celltype_NEP_condition")
VlnPlot(cancer_merge,features=c("dtd_P3NE35","dtd_MC7E37","dtd_P3CE36"),group.by = "celltype_NEP_condition")
RidgePlot(subset(cancer_merge,subset=celltype_NEP_condition=="PC3_40V_CytoD",invert=TRUE),features=c("dtd_FLD004","dtd_FLD500"),group.by = "celltype_NEP_condition")

FeaturePlot(cancer_merge,reduction="pca",
            features=c("CYBA","S100A2","PLAU","FOSL1"),max.cutoff = 3)
# FeaturePlot(subset(cancer_merge,subset=run=="first"),reduction="pca",
#               features=c("CYBA","S100A2","PLAU","FOSL1"),max.cutoff = 3)
# FeaturePlot(subset(cancer_merge,subset=run=="second"),reduction="pca",
#             features=c("CYBA","S100A2","PLAU","FOSL1"),max.cutoff = 3)
