library(ggrepel)
library(pheatmap)
library(destiny)
#https://satijalab.org/seurat/articles/integration_introduction.html

cancer_pc3 <- subset(subset(subset(cancer_merge123456,subset=celltype=="PC3"), subset=condition=="normal"),subset=NEP=="75V", invert = TRUE)
cancer_mda <- subset(cancer_merge123456,subset=celltype=="MDAMB231")
cancer_mcf7 <- subset(cancer_merge123456,subset=celltype=="MCF7")
cancer_mcf10a <- subset(subset(subset(cancer_merge123456,subset=celltype=="MCF10A"), subset=condition=="normal"), subset=run=="fifth", invert = TRUE)
#
#PC3 integrated
#
cancer <- cancer_pc3
cancer.list <- SplitObject(cancer, split.by = "NEP")
cancer.list <- lapply(X = cancer.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = cancer.list)
cancer.anchors <- FindIntegrationAnchors(object.list = cancer.list, anchor.features = features)
cancer.combined <- IntegrateData(anchorset = cancer.anchors)
DefaultAssay(cancer.combined) <- "integrated"
cancer.combined <- ScaleData(cancer.combined, verbose = FALSE)
cancer.combined <- RunPCA(cancer.combined, npcs = 30, verbose = FALSE)
cancer.combined <- RunUMAP(cancer.combined, reduction = "pca", dims = 1:30)
DimPlot(cancer.combined,reduction="umap" ,group.by = "NEP")+
  FeaturePlot(cancer.combined, features = "dtd_FLD004",reduction = "umap", max.cutoff = 4)
cancer_pc3.integrated <- cancer
rm(cancer, cancer_pc3)
#
#MDAMB231 integrated
#
cancer <- cancer_mda
cancer.list <- SplitObject(cancer, split.by = "NEP")
cancer.list <- lapply(X = cancer.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = cancer.list)
cancer.anchors <- FindIntegrationAnchors(object.list = cancer.list, anchor.features = features)
cancer.combined <- IntegrateData(anchorset = cancer.anchors)
DefaultAssay(cancer.combined) <- "integrated"
cancer.combined <- ScaleData(cancer.combined, verbose = FALSE)
cancer.combined <- RunPCA(cancer.combined, npcs = 30, verbose = FALSE)
cancer.combined <- RunUMAP(cancer.combined, reduction = "pca", dims = 1:30)
DimPlot(cancer.combined,reduction="umap" ,group.by = "NEP")+
  FeaturePlot(cancer.combined, features = "dtd_FLD004",reduction = "umap", max.cutoff = 4)
cancer_mda.integrated <- cancer
rm(cancer, cancer_mda)
#
#MCF7 integrated
#
cancer <- cancer_mcf7
cancer.list <- SplitObject(cancer, split.by = "run")
cancer.list <- lapply(X = cancer.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = cancer.list)
cancer.anchors <- FindIntegrationAnchors(object.list = cancer.list, anchor.features = features)
cancer.combined <- IntegrateData(anchorset = cancer.anchors)
DefaultAssay(cancer.combined) <- "integrated"
cancer.combined <- ScaleData(cancer.combined, verbose = FALSE)
cancer.combined <- RunPCA(cancer.combined, npcs = 30, verbose = FALSE)
cancer.combined <- RunUMAP(cancer.combined, reduction = "pca", dims = 1:30)
DimPlot(cancer.combined,reduction="umap" ,group.by = "NEP")+
  DimPlot(cancer.combined,reduction="umap" ,group.by = "run")+
  FeaturePlot(cancer.combined, features = "dtd_FLD004",reduction = "umap", max.cutoff = 4)
VlnPlot(cancer.combined, features = "dtd_FLD004", group.by = "run")+
  VlnPlot(cancer.combined, features = "dtd_FLD004", group.by = "NEP")
cancer_mcf7.integrated <- cancer
rm(cancer, cancer_mcf7)
#
#MCF10A integrated
#
cancer <- cancer_mcf10a
cancer.list <- SplitObject(cancer, split.by = "run")
cancer.list <- lapply(X = cancer.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = cancer.list)
cancer.anchors <- FindIntegrationAnchors(object.list = cancer.list, anchor.features = features)
cancer.combined <- IntegrateData(anchorset = cancer.anchors)
DefaultAssay(cancer.combined) <- "integrated"
cancer.combined <- ScaleData(cancer.combined, verbose = FALSE)
cancer.combined <- RunPCA(cancer.combined, npcs = 30, verbose = FALSE)
cancer.combined <- RunUMAP(cancer.combined, reduction = "pca", dims = 1:30)
DimPlot(cancer.combined,reduction="umap" ,group.by = "NEP")+
  DimPlot(cancer.combined,reduction="umap" ,group.by = "run")+
  FeaturePlot(cancer.combined, features = "dtd_FLD004",reduction = "umap", max.cutoff = 4)
cancer_mcf10a.integrated <- cancer
rm(cancer, cancer_mcf10a)
#
#Combine
#
cancer.integrated1 <-merge(cancer_pc3.integrated, y=cancer_mda.integrated)
cancer.integrated2 <-merge(cancer_mcf7.integrated, y=cancer_mcf10a.integrated)
cancer.integrated <-merge(cancer.integrated1, y=cancer.integrated2)
rm(cancer.list, cancer.combined, cancer.anchors, cancer.integrated1, cancer.integrated2)

cancer.integrated <- FindVariableFeatures(cancer.integrated, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer.integrated)
cancer.integrated <- ScaleData(cancer.integrated, features = all.genes)
cancer.integrated <- RunPCA(cancer.integrated, npcs=20, features = VariableFeatures(object = cancer.integrated))
DimPlot(cancer.integrated,group.by = "celltype")
cancer.integrated <- RunUMAP(cancer.integrated, dims = 1:10)
cancer.integrated <- RunTSNE(cancer.integrated,dims=1:10)
rm(cancer_pc3.integrated, cancer_mda.integrated, cancer_mcf7.integrated, cancer_mcf10a.integrated)


