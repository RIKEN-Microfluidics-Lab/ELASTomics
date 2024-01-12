### Required library ###
require(readr)
library(plyr)
library(dplyr,tidyr)
library(ggplot2)
library(tidyverse,R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(Seurat)
library(seqinr)
library(stringr)
library(ggrepel)
library(pheatmap)
library(destiny)
library(glmnet)
library("org.Hs.eg.db")
library(clusterProfiler)

### load ELASTomics data ###
load.elast.data <- function(wdir,cellname,min.cell.num){
  data <- Read10X(data.dir = wdir)
  tig <- CreateSeuratObject(data$`Gene Expression`,min.cells = min.cell.num)
  if (!is_empty(data$Custom)){
    tig[["DTD"]]<- CreateAssayObject(counts = data$Custom)
  }
  if (!is_empty(data$`Antibody Capture`)){
    tig[["ADT"]]<- CreateAssayObject(counts = data$`Antibody Capture`)
  }
  cellids <-colnames(tig)
  cellids<-paste0(cellname,substr(cellids,1,16))
  tig<-RenameCells(tig, new.names = cellids)
  return(tig)
}

# None-electroporated cell data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/CNTRL/outs/filtered_feature_bc_matrix/"
ctl <-load.elast.data(wdir,"CTL-",100)
ctl[["condition"]]<-"CTL"

# Electroporated cell  data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/EP/outs/filtered_feature_bc_matrix/"
nep <-load.elast.data(wdir,"NEP-",100)
nep[["condition"]]<-"NEP"

# Merge data
tig <- merge(ctl, y=nep)
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)


### Remove dead cells and multiple cells data ###
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^MT-")
tig <- subset(tig, subset= percent.mt<5)
tig <- subset(tig, subset = nCount_RNA > 2000)
tig <- subset(tig, subset = nCount_RNA < 60000)
FeatureScatter(tig,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA", group.by = "condition")


### PCA and visualize ###
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 300)
tig <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca",group.by = "condition")


### Identify the cell type ###
tig <- JackStraw(tig, num.replicate = 100)
tig <- ScoreJackStraw(tig, dims = 1:20)
JackStrawPlot(tig, dims = 1:20)
ElbowPlot(tig)
tig <- FindNeighbors(tig, dims = 1:19)
tig <- FindClusters(tig, resolution = 0.1)
cluster = as.numeric(Idents(tig))
tig <- RunUMAP(tig, dims = 1:19)
DimPlot(tig, reduction = "umap",group.by = "condition")+
  DimPlot(tig, reduction = "umap")+
  FeaturePlot(tig,features = c("ITGA8","DLGAP5"), reduction = "umap")
# DLGAP5: tig-1-20 (Young TIG-1 cells); ITGA8: tig-1-50 (Senescent TIG-1 cells)
tig <- RenameIdents(tig, `0` = "TIG1-50", `1` = "TIG1-20", `2` = "Unknown")
tig <- subset(tig,idents = "Unknown",invert=TRUE)

### Batch processing ###
#https://satijalab.org/seurat/articles/integration_introduction.html
tig.list <- SplitObject(tig, split.by = "condition")
tig.list <- lapply(X = tig.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = tig.list)
tig.anchors <- FindIntegrationAnchors(object.list = tig.list, anchor.features = features)
tig.combined <- IntegrateData(anchorset = tig.anchors)
DefaultAssay(tig.combined) <- "integrated"
tig.combined <- ScaleData(tig.combined, verbose = FALSE)
tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindClusters(tig.combined, resolution = 0.1)
tig.combined <- RenameIdents(tig.combined, `0` = "TIG1-50", `1` = "TIG1-20")
tig.combined[["age"]] <- Idents(tig.combined)
tig.combined <- NormalizeData(tig.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig.combined <- NormalizeData(tig.combined, normalization.method = "LogNormalize", scale.factor = 1e5)
tig.combined.nep <- subset(tig.combined,subset=condition=="NEP")
tig.combined.ctl <- subset(tig.combined,subset=condition=="CTL")
tig.combined <- RunUMAP(tig.combined, dims = 1:10)


### Output the analyzed images ###
# Fig. 3b
VlnPlot(tig.combined.nep,features = "dtd_FLD004", pt.size = 0.1)+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
# Fig. 3c
DimPlot(tig.combined, reduction = "umap", label = FALSE,group.by = "condition", cols = c("#0072B2","#D55E00"))
# Fig. 3d
DimPlot(tig.combined, reduction = "umap", label = FALSE, cols = c("#0072B2","#D55E00"))
# Fig. 3e
FeaturePlot(tig.combined.nep, reduction = "umap", feature = 'FLD004', slot = "data")
# Supplementary Fig. 12a
VlnPlot(tig.combined.ctl,features = "CDKN1A")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "DCN")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "TIMP1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "LMNB1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "HMGB1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "CAV1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
# Supplementary Fig. 12b
VlnPlot(tig.combined.nep,features = "CDKN1A")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "DCN")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "TIMP1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "LMNB1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "HMGB1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "CAV1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
# Supplementary Fig. 12c
FeatureScatter(tig.combined.nep, feature1 = 'FLD004', feature2 = 'FOS', cols = c("#0072B2","#D55E00"))
# Supplementary Fig. 12d
FeatureScatter(tig.combined.nep, feature1 = 'FLD004', feature2 = 'RRAD', cols = c("#0072B2","#D55E00"))
# Supplementary Fig. 12e
FeatureScatter(tig.combined.nep, feature1 = 'FLD004', feature2 = 'ABL2', cols = c("#0072B2","#D55E00"))
# Supplementary Fig. 12f
FeatureScatter(tig.combined.nep, feature1 = 'FLD004', feature2 = 'YWHAH', cols = c("#0072B2","#D55E00"))


### Correlation ###
exp.matrix <- data.frame(t(data.frame(tig.combined.nep[["integrated"]]@data)))
var.gene <- data.frame(VariableFeatures(tig.combined.nep))
response <- data.frame(t(data.frame(tig.combined.nep[["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
cors$bool <- FALSE
cors[cors$gene %in% c("FOS", "ADM", "CDKN1A", "DCN", "CITED2", "TIMP1", "PTGS2", "NQO1", "IGFBP5", "CCN2", "SOD2", "CTSC", "PLK2", "EDNRB"),]$bool<-TRUE
cors[cors$gene %in% c("DLGAP5", "NCL", "PCLAF", "DTYMK"),]$bool<-TRUE
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.15)
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
# Fig. 3f
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.15),aes(x=rank,y=cors, color = "red"))+theme_classic() +NoLegend()+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1) 


### Glmnet ###
#https://glmnet.stanford.edu/articles/glmnet.html
#https://www.r-bloggers.com/2021/05/lasso-regression-model-with-r-code/
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
exp.matrix <- t(data.frame(tig.combined.nep[["integrated"]]@data))
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
alpha <- seq(0, 1, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = best.alpha)
best.lambda <- m$lambda.min
en.model <- glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian",
                   lambda = best.lambda, alpha = best.alpha)
en.model.beta <-data.frame(coef(en.model, s="lambda.min"))
en.model.beta$gene <- rownames(en.model.beta)
en.model.beta$cors <- cors.sig[en.model.beta$gene,]$cors
en.model.beta <- na.omit(en.model.beta)
en.model.plus.beta <- en.model.beta[en.model.beta$s1>0,]
en.model.minus.beta <- en.model.beta[en.model.beta$s1<0,]
# Fig. 3g
ggplot(en.model.beta,aes(x=cors,y=s1))+geom_point()+
  geom_point(data=subset(en.model.beta,gene %in% c("AC007952.4","AC091271.1","AC021155.5","KLF2","FOS","RRAD", "KLF4", "IGFBP4", "ADM", "BTG2", "RPL22L1","YWHAH", "SMC3", "ABL2")),aes(y=s1,x=cors,color="red"))+
  geom_text_repel(data=subset(en.model.beta,gene %in% c("AC007952.4","AC091271.1","AC021155.5","KLF2","FOS","RRAD", "KLF4", "IGFBP4", "ADM", "BTG2", "RPL22L1","YWHAH", "SMC3", "ABL2")),aes(y=s1,x=cors,label=gene,fontface = "italic"))+
  theme_bw()+NoLegend()

### GSEA ###
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[gene.list$symbol,1] #correlation
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
gene_list_log2fc_cor <- symbol2entrez.order(cors)
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       verbose      = FALSE)
View(gse_result_cor@result)

# Supplementary Fig. 12g
ridgeplot(gse_result_cor,showCategory = 40)
# Supplementary Fig. 12h-k
g1 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0140657") #ATP-dependent activity
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0016887") #ATP hydrolysis activity
g3 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0003779") #actin binding
g4 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0007568") #aging
gridExtra::grid.arrange(g3, g4, nrow = 2) 
