library(stringr)
library(Seurat)
library(ggplot2)
library("org.Hs.eg.db")
#datadir <-"/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_TIGbulk/"
#rdir <- "/home/samba/public/shintaku/github/ELASTomics/"
wdir <-"/home/samba/public/shintaku/ELASTomics/20220719HiSeqX014_TIGbulk/"
barcode <- read.table(file.path("/home/samba/crick/shintaku/github/hunter2/cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

#source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'unused/preprocess/preprocess_whitelist.R'))
active_barcode <- barcode#barcode[sort(unique(allencoded$first_index)),]
encode_barcode=TRUE
source("unused/preprocess/preprocess_RNAseq_data.R")
source("unused/preprocess/preprocess_save_10x_format.R")


tig.bulk <- Read10X(data.dir = wdir)

annotate.samples <-function(barcodes,PDL){
  cell <- barcodes
  cell$celltype <- PDL
#  cell$culture <- culture
  colnames(cell)<-c("RTid","PDL")
  return(cell)
}


tig.pdl.33<-annotate.samples(data.frame(c("1","5","7")),"PDL33")
tig.pdl.44<-annotate.samples(data.frame(c("9","14","17")),"PDL44")
tig.pdl.52 <- annotate.samples(data.frame(c("19","21","24")),"PDL52")
tig.pdl.57<-annotate.samples(data.frame(c("28","29","30")),"PDL57")
tig.pdl.65<-annotate.samples(data.frame(c("31","33","36")),"PDL65")
sample.meta.data <- rbind(tig.pdl.33,tig.pdl.44,tig.pdl.52,tig.pdl.57,tig.pdl.65)
rownames(sample.meta.data)<-sample.meta.data$RTid

cellids <- colnames(tig.bulk)#substr(colnames(tig.bulk),12,22)
#colnames(tig.bulk)<-cellids

tig.bulk.seurat <- CreateSeuratObject(tig.bulk)
tig.bulk.seurat[["batch"]]<-substr(cellids,9,10)
tig.bulk.seurat[["RTid"]]<-substr(cellids,12,13)
tig.bulk.seurat[["PDL"]]<-sample.meta.data[substr(colnames(tig.bulk.seurat),12,13),]$PDL
#tig.bulk.seurat[["culture"]]<-sample.meta.data[colnames(tig.bulk.seurat),]$culture

FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "PDL")+
  FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "batch")

tig.bulk.seurat <- NormalizeData(tig.bulk.seurat, normalization.method = "LogNormalize", scale.factor = 1e5)
tig.bulk.seurat <- FindVariableFeatures(tig.bulk.seurat, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(tig.bulk.seurat)

tig.bulk.seurat <- ScaleData(tig.bulk.seurat, features = all.genes)
tig.bulk.seurat <- RunPCA(tig.bulk.seurat, npcs=3, features = VariableFeatures(object = tig.bulk.seurat))
DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "PDL")+
  DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "batch")


#tig.bulk.seurat[["type_culture"]]<- paste0(unlist(tig.bulk.seurat[["celltype"]]),"_",unlist(tig.bulk.seurat[["culture"]]))
tig.bulk.seurat <- subset(tig.bulk.seurat,subset=batch=="03")
#Idents(tig.bulk.seurat)<-tig.bulk.seurat[["culture"]]
Idents(tig.bulk.seurat)<-tig.bulk.seurat[["PDL"]]

#tig.marker <- FindMarkers(subset(tig.bulk.seurat,subset=culture=="mono"),ident.1 = "tig50",ident.2 = "tig20", test.use = "wilcox")
#tig.marker$gene <-rownames(tig.marker)
#View(tig.marker)
#https://encycolorpedia.com/d55e00
#https://encycolorpedia.jp/0072b2
p1<-VlnPlot(tig.bulk.seurat, features = "ATF3",group.by = "PDL")+ylim(c(0,5))
p2<-VlnPlot(tig.bulk.seurat, features = "FOS",group.by = "PDL")+ylim(c(0,5))
p3<-VlnPlot(tig.bulk.seurat, features = "ADM",group.by = "PDL")+ylim(c(0,5))
p4<-VlnPlot(tig.bulk.seurat, features = "CDKN1A",group.by = "PDL")+ylim(c(0,5))
p6<-VlnPlot(tig.bulk.seurat, features = "DTYMK",group.by = "PDL")+ylim(c(0,5))
p7<-VlnPlot(tig.bulk.seurat, features = "NCL",group.by = "PDL")+ylim(c(0,5))
p8<-VlnPlot(tig.bulk.seurat, features = "PCLAF",group.by = "PDL")+ylim(c(0,5))
p9<-VlnPlot(tig.bulk.seurat, features = "DLGAP5",group.by = "PDL")+ylim(c(0,5))
p5<-VlnPlot(tig.bulk.seurat, features = "RRAD",group.by = "PDL")+ylim(c(0,5))
p10<-VlnPlot(tig.bulk.seurat, features = "ACTB",group.by = "PDL")+ylim(c(0,6))
gridExtra::grid.arrange(p1, p2,p3, p4, p5, p6,p7,p8,p9,p10,nrow = 2)


p1<-VlnPlot(tig.bulk.seurat, features = "RRAD",group.by = "PDL")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#D55E00","#ED945E"))
p2<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig50"), features = "RRAD",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#0072B2","#769ECC"))
p3<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig20"), features = "RHOB",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#D55E00","#ED945E"))
p4<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig50"), features = "RHOB",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#0072B2","#769ECC"))
p5<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig20"), features = "ACTB",group.by = "culture")+ylim(c(0,5))+
  scale_fill_manual(values = c("#D55E00","#ED945E"))
p6<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig50"), features = "ACTB",group.by = "culture")+ylim(c(0,5))+
  scale_fill_manual(values = c("#0072B2","#769ECC"))
p7<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig20"), features = "GADD45B",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#D55E00","#ED945E"))
p8<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig50"), features = "GADD45B",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#0072B2","#769ECC"))
p9<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig20"), features = "ABL2",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#D55E00","#ED945E"))
p10<-VlnPlot(subset(tig.bulk.seurat,subset=celltype=="tig50"), features = "ABL2",group.by = "culture")+ylim(c(0,2.5))+
  scale_fill_manual(values = c("#0072B2","#769ECC"))
gridExtra::grid.arrange(p1, p2,p3, p4, p5, p6,p7,p8,p9,p10,nrow = 5)

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
#
#  plot senescence levels of TIG-1 cell populations
#
age <- tig.bulk.seurat[["celltype"]]
culture <- tig.bulk.seurat[["culture"]]
#colnames(age)<-c("celltype","culture")
cellids<-colnames(tig.bulk.seurat)
age<-data.frame(paste0(unlist(age),".",unlist(culture)))
rownames(age)<-colnames(tig.bulk.seurat)
colnames(age)<-"sample"
age$celltype<-tig.bulk.seurat[["celltype"]]

exp.matrix<-data.frame(tig.bulk.seurat[["RNA"]]@data)
exp.matrix.age<-exp.matrix[rownames(exp.matrix) %in% aging_gene,]
#cors.sig.sub <-cors.sig[!(cors.sig$gene %in% "RRAD"),]

exp.matrix.cor<-exp.matrix[rownames(exp.matrix) %in% cors.sig$gene,]
dm.age <- DiffusionMap(data.frame(t(exp.matrix.age)))
dm.cor <- DiffusionMap(data.frame(t(exp.matrix.cor)))
age$dm.age<-dm.age$DC1
age$dm.cor<-dm.cor$DC1

p1<-ggplot(age,aes(x=factor(sample,levels=c("tig20.mono","tig20.co","tig50.co","tig50.mono")),y=dm.age,fill=sample))+
  geom_violin()+geom_jitter()+
  scale_fill_manual(values = c("#D55E00","#ED945E", "#0072B2", "#769ECC"))+ theme_bw()
p2<-ggplot(age,aes(x=factor(sample,levels=c("tig20.mono","tig20.co","tig50.co","tig50.mono")),y=dm.cor,fill=sample))+
  geom_violin()+geom_jitter()+
  scale_fill_manual(values = c("#D55E00","#ED945E", "#0072B2", "#769ECC"))+ theme_bw()

gridExtra::grid.arrange(p1, p2, nrow = 2)


pheatmap(exp.matrix[,order(dm$DC1,decreasing = FALSE)],
         scale="row",
         cluster_cols = FALSE,
         annotation_col = age,
         show_colnames = FALSE)
age$DC1<-dm$DC1
ggplot(age,aes(x=age,y=-DC1))+geom_violin()+geom_jitter(size=0.1)
