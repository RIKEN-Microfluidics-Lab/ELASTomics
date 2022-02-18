library(stringr)
library(Seurat)
library(ggplot2)
tig.bulk <- Read10X(data.dir = "/home/samba/sanger/shintaku/20220109HiSeqX008_tig/")

annotate.samples <-function(barcodes,celltype,culture){
  cell <- barcodes
  cell$celltype <- celltype
  cell$culture <- culture
  colnames(cell)<-c("barcode","celltype","culture")
  return(cell)
}
tig.20.mono<-annotate.samples(data.frame(c("ACCTCGTTGA","CCATTGCGTT","ATCTCCGAAC")),"tig20","mono")
tig.50.mono<-annotate.samples(data.frame(c("GGCAGGTATT","AATCGTAGCG","GGTACTGCCT")),"tig50","mono")
tig.20.co <- annotate.samples(data.frame(c("GCCATTCTCC","TGGTCAGCCA","TCGGCCTTAC")),"tig20","co")
tig.50.co<-annotate.samples(data.frame(c("TTACCGAGGC","GATACGGAAC","AGAACGTCTC")),"tig50","co")
sample.meta.data <- rbind(tig.20.mono,tig.50.mono,tig.20.co,tig.50.co)
rownames(sample.meta.data)<-sample.meta.data$barcode

cellids <- substr(colnames(tig.bulk),12,22)
colnames(tig.bulk)<-cellids

tig.bulk.seurat <- CreateSeuratObject(tig.bulk)
tig.bulk.seurat[["celltype"]]<-sample.meta.data[colnames(tig.bulk.seurat),]$celltype
tig.bulk.seurat[["culture"]]<-sample.meta.data[colnames(tig.bulk.seurat),]$culture

FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "celltype")+
  FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "culture")

tig.bulk.seurat <- NormalizeData(tig.bulk.seurat, normalization.method = "LogNormalize", scale.factor = 1e5)
tig.bulk.seurat <- FindVariableFeatures(tig.bulk.seurat, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(tig.bulk.seurat)

tig.bulk.seurat <- ScaleData(tig.bulk.seurat, features = all.genes)
tig.bulk.seurat <- RunPCA(tig.bulk.seurat, npcs=3, features = VariableFeatures(object = tig.bulk.seurat))
DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "celltype")+
DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "culture")

Idents(object = tig.bulk.seurat) <- tig.bulk.seurat[["celltype"]]

tig.bulk.seurat<- CellCycleScoring(tig.bulk.seurat,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
DimPlot(tig.bulk.seurat,reduction="pca",group.by = "celltype")+
  DimPlot(tig.bulk.seurat,reduction="pca")


tig.mono <-subset(tig.bulk.seurat,subset=culture=="mono")
tig.de.markers <- FindMarkers(tig.bulk.seurat, ident.1 = "tig20", ident.2 = "tig50")
tig.de.markers[ abs(tig.de.markers$avg_log2FC)>1.0,]

ggplot(tig.de.markers,aes(y=-log(p_val),x=avg_log2FC))+geom_point()

DimPlot(tig.bulk.seurat, reduction="pca", group.by = "celltype") + 
  FeaturePlot(tig.bulk.seurat,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"))
