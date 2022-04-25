library(stringr)
library(Seurat)
library(ggplot2)
tig.bulk <- Read10X(data.dir = "/home/samba/sanger/shintaku/ELASTomics/20220109HiSeqX008_tig/")

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


tig.bulk.seurat[["type_culture"]]<- paste0(unlist(tig.bulk.seurat[["celltype"]]),"_",unlist(tig.bulk.seurat[["culture"]]))

Idents(tig.bulk.seurat)<-tig.bulk.seurat[["type_culture"]]

tig.marker <- FindMarkers(tig.combined.ctl,ident.1 = "TIG1-50",ident.2 = "TIG1-20", test.use = "wilcox",
                          print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf)

symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}


gene.list<-symbol2entrez(en.model.minus.beta$gene)

gene.list<-symbol2entrez(rownames(tig.marker[tig.marker$avg_log2FC>0.2 & tig.marker$p_val_adj<0.1,]))
ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))


library("org.Hs.eg.db")
library(clusterProfiler)
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  #gene.list$s0 <- genes[gene.list$symbol,1]
  gene.list$s0 <- genes[gene.list$symbol,1]
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}

genes<-data.frame(tig.marker$avg_log2FC)
rownames(genes)<-rownames(tig.marker)
gene_list_log2fc <- symbol2entrez.order(genes)
#
# try BP: biological process, CC: cellular component, or MF: molecular function
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 12,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
ridgeplot(gse_result,showCategory = 20)
View(gse_result@result)

gene_list <- unlist(strsplit(gse_result@result[gse_result@result$Description=="sister chromatid segregation",]$core_enrichment,"/"))
SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
gene_list_name <- SYBOL2EG[SYBOL2EG$gene_id %in% gene_list,]




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
