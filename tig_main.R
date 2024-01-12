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
rdir <- "/home/samba/public/shintaku/github/ELASTomics/"

# load ELASTomics data from an output of cellranger with cite-seq pipeline
source("elast.load.elast.data.R")

# control data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/CNTRL/outs/filtered_feature_bc_matrix/"
ctl <-load.elast.data(wdir,"CTL-",100)
ctl[["condition"]]<-"CTL"

# elastomics data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/EP/outs/filtered_feature_bc_matrix/"
nep <-load.elast.data(wdir,"NEP-",100)
nep[["condition"]]<-"NEP"

tig <- merge(ctl, y=nep)
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)

# mitochondrial gene percent and remove dead cells
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^MT-")
tig <- subset(tig, subset= percent.mt<5)
tig <- subset(tig, subset = nCount_RNA > 2000)
tig <- subset(tig, subset = nCount_RNA < 60000)
FeatureScatter(tig,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA", group.by = "condition")

# merge all the tig data
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 300)
top10 <- head(VariableFeatures(tig), 40)
plot1 <- VariableFeaturePlot(tig)
plot1 <- LabelPoints(plot = tig, points = top10)
plot1
# normalize the dtd data with "RC" option.
tig <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)

# PCA and visualize
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca",group.by = "condition")

source("Seurat.clustering.R")

#subset elastomics data after the normalization
tig.nep<-subset(tig,subset=condition=="NEP")
tig.ctl<-subset(tig,subset=condition=="CTL")
DimPlot(tig.nep, reduction = "pca", label = TRUE)

# integrate elast result into Seurat object
source("elast.integrated.R")
tig.combined.nep <- subset(tig.combined,subset=condition=="NEP")
tig.combined.ctl <- subset(tig.combined,subset=condition=="CTL")
tig.combined <- RunUMAP(tig.combined, dims = 1:10)
DimPlot(tig.combined, reduction = "umap", label = FALSE, cols = c("#0072B2","#D55E00"))
DimPlot(tig.combined, reduction = "umap", label = FALSE,group.by = "condition", cols = c("#0072B2","#D55E00"))
FeaturePlot(tig.combined.nep, reduction = "umap", feature = 'FLD004', slot = "data")
FeaturePlot(tig.combined.nep, reduction = "umap", feature = 'nFeature_RNA', slot = "data")
VlnPlot(tig.combined.nep,features = "dtd_FLD004", pt.size = 0.1)+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "CDKN1A")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c("#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.nep,features = "CAV1")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
VlnPlot(tig.combined.ctl,features = "CDKN1A")+scale_x_discrete(limits=c("TIG1-20", "TIG1-50"))+scale_fill_manual(values = c( "#0072B2","#D55E00"))+NoLegend()
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
FeaturePlot(tig.combined.nep, reduction = "pca", feature = 'FLD004', slot = "data")
FeatureScatter(tig.combined.nep, feature1 = 'FLD004', feature2 = 'YWHAH', cols = c("#0072B2","#D55E00"))
summary(tig.combined.nep$nCount_RNA)
summary(tig.combined.nep$nFeature_RNA)

source("elast.glmnet.R")
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  #  geom_text(aes(label=ifelse(cors>0.2,as.character(gene),'')),hjust=0, vjust=0)+
  geom_point(data=subset(cors, gene %in% c("AC007952.4","AC091271.1","ABL2","FOS","RRAD","RHOB","RPL22L1","YWHAH","RASD1","AL021155.5","BTG2","KLF2")),aes(y=cors,x=rank,color="red"))+
  geom_text_repel(data=subset(cors,gene %in% c("AC007952.4","AC091271.1","ABL2","FOS","RRAD","RHOB","RPL22L1","YWHAH","RASD1","AL021155.5","BTG2","KLF2")),aes(y=cors,x=rank,label=gene,fontface = "italic"))+
  theme_bw()+ylim(c(-0.25,0.5))+NoLegend()



#GseGo
MCF7 <- FindMarkers(tig.combined.nep,ident.1 = "TIG1-50",ident.2 = "TIG1-20", test.use = "wilcox", , min.pct = 0.25, logfc.threshold = 0)
MCF7$gene <- rownames(MCF7)
MCF7 <- inner_join(MCF7, cors, by="gene")
rownames(MCF7) <- MCF7$gene 
source("cancer.GSEA.R")

#cytoskeletal protein binding
g1 <- gseaplot(gse_result_fol, by = "runningScore", geneSetID = "GO:0008092")
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0008092")
gridExtra::grid.arrange(g1, g2, nrow = 2) 
#ATP hydrolysis activity
g1 <- gseaplot(gse_result_fol, by = "runningScore", geneSetID = "GO:0016887")
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0016887")

#Comparison
library(gplots)
Elec <- FindMarkers(tig.combined, group.by="condition", ident.1 = "NEP", ident.2 = "CTL", min.pct = 0.25, logfc.threshold = 0)
Elec$gene <- rownames(Elec)
TIG1 <- inner_join(Elec, cors, by="gene")
ggplot(TIG1,aes(x=avg_log2FC, y=cors))+geom_point(color = "black")+theme_classic()
ggplot(TIG1,aes(x=avg_log2FC, y=cors))+geom_point()+
  geom_point(data=subset(TIG1, abs(avg_log2FC) > 2.5),aes(x=avg_log2FC,y=cors, color = "red"))+xlim(-10, 10)+theme_classic()




# re-scale the dtd data with the concentration of the dtd molecules in the solution
source("elast.rescale.dtd.R")
# FLD004,FLD010,FLD040,FLD070,FLD150,FLD500, and others
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "RC",scale.factor = 1e2)
concentration<- data.frame(c(2,6,9,3,9,3,1,1,1,1,1,1,1))
tig.combined.nep<-rescale.dtd.data(tig.combined.nep,concentration)
tig.dtd.scale<-tig.combined.nep[["DTD"]]@data
#tig.dtd.scale <- tig.dtd.scale[c("FLD004","FLD010","FLD070","FLD500"),]
# compute a radius for a single case
source("elast.comp.radii.R")
#
# S.radii: input Stokes radii of DTD
# tig.dtd.scale; scaled amount of DTD imported to cells
# The nrow of the tig.dtd.scale must match with the length of S.radii
#
S.radii <- c(4.1e-9, 4.2e-9, 7.8e-9, 10.6e-9, 15.1e-9, 17.0e-9) # Stokes radii of DTD
#S.radii <- c(1.4e-9, 2.4e-9, 1.4e-9, 1.4e-9, 1.4e-9, 1.4e-9)
#fm <-elast.comp.radius(tig.dtd.scale[1:4,"NEP-CCCTTAGGTCAAACGG"],S.radii,TRUE)
#
# compute radii for multiple cases
res<-elast.comp.radii(tig.dtd.scale[1:4,],S.radii,TRUE)
Res <- res[[1]]
dtd.predict<-res[[2]]
#
# extract explanatory variables via glmnet
#





en.model.nonzero.beta <- subset(en.model.beta, subset=abs(s0)>0.001)
#source("elast.biomaRt.R")
en.model.plus.beta <- subset(en.model.beta, subset=s0> 0.001)
en.model.minus.beta <- subset(en.model.beta, subset=s0< -0.001)



source(file.path(rdir,"elast.integrated.cor.diffusionmap.R"))

#GSEA
exp.matrix <- t(data.frame(tig.combined.nep [["RNA"]]@data))
var.gene <- data.frame(VariableFeatures(tig.combined.nep ))
response <- data.frame(t(data.frame(tig.combined.nep [["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
cors$bool <- FALSE
cors[cors$gene %in% rownames(cors %>% top_n(5, cors)),]$bool<-TRUE
cors[cors$gene %in% rownames(cors %>% top_n(-5, cors)),]$bool<-TRUE
cors.sig<-subset(subset(cors,subset=pval<1e-4),subset=abs(cors)>0.2)
cors.pos<-subset(subset(cors,subset=pval<1e-4),subset = cors > 0.2)
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.2),aes(x=rank,y=cors, color = "red"))+theme_classic() +ylim(-0.35, 0.35) +NoLegend()+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1) 
library("org.Hs.eg.db")
library(clusterProfiler)
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
ridgeplot(gse_result_cor,showCategory = 40)

g1 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0140657") #ATP-dependent activity
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0016887") #ATP hydrolysis activity
g3 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0003779") #actin binding
g4 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0007568") #aging
gridExtra::grid.arrange(g3, g4, nrow = 2) 
