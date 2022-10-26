require(readr)
library(plyr)
library(dplyr,tidyr)
library(ggplot2)
library(tidyverse,R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(Seurat)
library(SingleCellSignalR)
library(seqinr)
library(stringr)
rdir <- "/home/samba/public/shintaku/github/ELASTomics/"
#
# load ELASTomics data from an output of cellranger with cite-seq pipeline
#
source("elast.load.elast.data.R")
# control data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_10x/mhsc_1_14c_2/outs/filtered_feature_bc_matrix/"
tig <-load.elast.data(wdir,"NEP-",100)
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)
# mitochondrial gene percent and remove dead cells
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^mt-")
tig <- subset(tig, subset= percent.mt<5)
tig <- subset(tig,subset=nCount_RNA >1000)
#
#
# merge all the tig data
#
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 300)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tig), 40)
# plot variable features with labels
plot1 <- VariableFeaturePlot(tig)
plot1 <- LabelPoints(plot = tig, points = top10)
plot1

#reject Abnormal replication


# normalize the dtd data with "RC" option.
#tig <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig <- NormalizeData(tig,assay="ADT",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca")
source("Seurat.clustering.R")
JackStrawPlot(tig, dims = 1:20)
ElbowPlot(tig)

#Celltype
tig <- FindNeighbors(tig, dims = 1:18)
tig <- FindClusters(tig, resolution = 0.8)
DimPlot(tig, label = FALSE)
new.cluster.ids <- c("Gra", "Gra", "Gra","Gra", "GMP", "Mono","Gra", "EryP", "MMP","Erythroid", "CMP", "Bcell","Bcell", "Erythroid", "CLP","Gra", "CMP", "Other","Other")
names(new.cluster.ids) <- levels(tig)
tig <- RenameIdents(tig, new.cluster.ids)
tig[["Celltype"]] <- Idents(tig)

#CMP or MkP
tig1 <- subset(tig, idents = c("CMP"))
tig1 <- FindNeighbors(tig1, dims = 1:10)
tig1 <- FindClusters(tig1, resolution = 0.2)
DimPlot(tig1, label = FALSE)
new.cluster.ids <- c("MkP", "CMP", "CMP","CMP")
names(new.cluster.ids) <- levels(tig1)
tig1 <- RenameIdents(tig1, new.cluster.ids)
tig1[["Celltype"]] <- Idents(tig1)

#HSC or MMP
tig2 <- subset(tig, idents = c("MMP"))
tig2 <- FindNeighbors(tig2, dims = 1:20)
tig2 <- FindClusters(tig2, resolution = 0.8)
DimPlot(tig2, label = FALSE)
new.cluster.ids <- c("HSC", "MMP", "MMP", "MMP", "MMP")
names(new.cluster.ids) <- levels(tig2)
tig2 <- RenameIdents(tig2, new.cluster.ids)
tig2[["Celltype"]] <- Idents(tig2)

#Reconstruction
tig3 <- subset(tig, idents = c("CMP", "MMP"), invert = TRUE)
tig12 <- merge(tig1, y=tig2)
tig4 <- merge(tig3, y=tig12)
tig4 <- FindVariableFeatures(tig4, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(tig4)
tig4 <- ScaleData(tig4, features = all.genes)
tig4 <- RunPCA(tig4, npcs=20, features = VariableFeatures(object = tig4))
#tig4 <- FindNeighbors(tig4, dims = 1:19)
#tig4 <- FindClusters(tig4, resolution = 0.1)
tig4 <- RunUMAP(tig4, dims = 1:19)
DimPlot(tig, reduction = "umap", label =TRUE, pt.size = 0.5, group.by = "Celltype") #+ NoLegend()
tig <- tig4
rm(tig1, tig2, tig12, tig3, tig4)

#Gene expression
DimPlot(tig, reduction = "umap", label =TRUE, pt.size = 0.5) #+ NoLegend()
FeaturePlot(tig,features = c("Camp", "Chil3","Mmp9"),reduction="umap", ncol = 3) #Granulocyte
FeaturePlot(tig,features = c("Elane", "Vcam1","Mpo"),reduction="umap", ncol = 3) #GMP
FeaturePlot(tig,features = c("Ccr2", "Klf4","Irf8"),reduction="umap", ncol = 3) #Mono
FeaturePlot(tig,features = c("Cox6a2", "Dntt","Il7r"),reduction="umap", ncol = 3) #CLP
FeaturePlot(tig,features = c("Ly6a", "Hlf", "Rgs1"),reduction="umap", ncol = 3) #MMP+HSC
VlnPlot(tig, features = c("Hlf"), log = TRUE,  slot = "count") #HSC
FeaturePlot(tig,features = c("Cd79a", "Vpreb1", "Pax5"),reduction="umap", ncol = 3) #other(B-cell)
FeaturePlot(tig,features = c("Itga2b", "Gata2", "Pf4"),reduction="umap", ncol = 3) #CMP+MkP
VlnPlot(tig, features = c("Pf4"), log = TRUE,  slot = "count") #MkP
FeaturePlot(tig,features = c("Klf1", "Mt2", "Gata1"),reduction="umap", ncol = 3) #EryP
FeaturePlot(tig,features = c("Gypa", "Slc4a1", "Hba-a1"),reduction="umap", ncol = 3) #Erythroid
FeaturePlot(tig,features = "percent.mt",reduction="umap")

#CITE-seq
tig <- subset(tig, subset= adt_ADT429<5000)
tig <- subset(tig, subset= adt_ADT130<5000)
FeaturePlot(tig,features="adt_ADT130",max.cutoff = 3) #ADT130 = Ly6a(Sca-1)
VlnPlot(tig, features = c("adt_ADT130"), log = TRUE,  slot = "count") 
FeaturePlot(tig,features="adt_ADT429",max.cutoff = 3) #ADT429 = CD48
VlnPlot(tig, features = c("adt_ADT429"), log = TRUE,  slot = "count")
FeaturePlot(tig,features="adt_ADT012",max.cutoff = 200) #ADT012 = c-kit(CD117)
VlnPlot(tig, features = c("adt_ADT012"), log = TRUE,  slot = "count")
FeaturePlot(tig,features="adt_ADT203",max.cutoff = 2.5) #ADT203 = CD150(SLAM)
VlnPlot(tig, features = c("adt_ADT203"), log = TRUE,  slot = "count")

#ELASTomics
FeaturePlot(tig,features="FLD150",max.cutoff = 10) 
VlnPlot(tig, features = c("dtd_FLD500"), log = TRUE,  slot = "counts")
FeatureScatter(object = tig, feature1 = 'adt_ADT130', feature2 = 'dtd_FLD500')
p <- pivot_longer(as.data.frame(t(tig[["DTD"]]@data)), cols = 1:6, names_to = "Categories", values_to = "Values")
g <- ggplot(p, aes(x = Categories, y = Values))
g+geom_violin()+scale_y_log10()+geom_jitter(width=0.3)


FeaturePlot(tig,features="Ly6a",max.cutoff = 2)
FeaturePlot(tig,features="dtd_FLD500",max.cutoff = 4)
VlnPlot(tig, features = c("dtd_FLD500"), log = TRUE,  slot = "data")
VlnPlot(tig, features = c("adt_ADT130"), log = TRUE,  slot = "count")




#Correlation
library(glmnet)
tig1 <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig1 <- NormalizeData(tig1,assay="ADT",normalization.method = "CLR",scale.factor = 1e2)
tig1 <- subset(tig1, subset= dtd_FLD500>1)
tig1 <- subset(tig1, subset= dtd_FLD500<5)
tig1 <- subset(tig1, idents = c("Other"), invert= TRUE)
RidgePlot(object = tig1, features='dtd_FLD500')
FeatureScatter(object = tig1, feature1 = 'adt_ADT130', feature2 = 'dtd_FLD500')
tig1 <- subset(tig1, idents = c("CMP", "EryP", "Erythroid"))


exp.matrix <- t(data.frame(tig1[["RNA"]]@data))
var.gene <- data.frame(VariableFeatures(tig1))
response <- data.frame(t(data.frame(tig1[["DTD"]]@data)))

rownames(response) <- rownames(exp.matrix)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD500,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
#cors$bool <- FALSE
#cors[cors$gene %in% "IGFBP2",]$bool<-TRUE
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()
FeatureScatter(object = tig1, feature1 = 'Spta1', feature2 = 'dtd_FLD500')
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.05)

#egoGo
library(clusterProfiler)
library(org.Mm.eg.db)
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}
#gene <- symbol2entrez(append(rownames(cors[cors$cors > 0.2,]),rownames(cors[cors$cors < -0.2,])))
gene <- symbol2entrez(rownames(cors.sig))
ego_result <- enrichGO(gene          = gene$gene_id, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))

# heatmap of correlated genes
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.05)
exp.matrix <- data.frame(tig1[["RNA"]]@data)
dtd <-FetchData(object = tig1, vars = c("dtd_FLD500"))
rownames(dtd)<-colnames(exp.matrix)
out<-pheatmap(exp.matrix[rownames(cors.sig),order(dtd$dtd_FLD500,decreasing = FALSE)],
              cluster_cols = FALSE, show_colnames = FALSE,
              annotation_col = dtd,
              clustering_distance_rows = "correlation")


# elastic net
exp.matrix <- t(data.frame(tig1[["RNA"]]@data))
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
alpha <- seq(0, 1, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = exp.matrix, y = response$FLD500, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
m <- cv.glmnet(x = exp.matrix, y = response$FLD500, family = "gaussian", alpha = best.alpha)
best.lambda <- m$lambda.min
en.model <- glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian",
                   lambda = best.lambda, alpha = best.alpha)
en.model.beta <-data.frame(coef(en.model, s="lambda.min"))
en.model.beta$gene <- rownames(en.model.beta)
en.model.beta$cors <- cors.sig[en.model.beta$gene,]$cors
en.model.beta <- na.omit(en.model.beta)
ggplot(en.model.beta, aes(x = cors, y = s1, label = gene))+geom_point()+theme_bw()+ylim(c(-0.25,0.25)) + #geom_text_repel()
  geom_text_repel(aes(label=ifelse(abs(s1) > 0.1,as.character(gene),'')))






#
# re-scale the dtd data with the concentration of the dtd molecules in the solution
# 
source("elast.rescale.dtd.R")
# FLD004,FLD010,FLD040,FLD070,FLD150,FLD500, and others
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "RC",scale.factor = 1e2)
concentration<- data.frame(c(2,6,9,3,9,3,1,1,1,1,1,1,1))
tig.combined.nep<-rescale.dtd.data(tig.combined.nep,concentration)

tig.dtd.scale<-tig.combined.nep[["DTD"]]@data
#tig.dtd.scale <- tig.dtd.scale[c("FLD004","FLD010","FLD070","FLD500"),]


#
# compute a radius for a single case
source("elast.comp.radii.R")
#
# S.radii: input Stokes radii of DTD
# tig.dtd.scale; scaled amount of DTD imported to cells
# The nrow of the tig.dtd.scale must match with the length of S.radii
#
S.radii <- c(1.4e-9, 2.3e-9, 4.5e-9, 6.3e-9, 8.5e-9, 15.1e-9) # Stokes radii of DTD

#S.radii <- c(1.4e-9, 2.4e-9, 1.4e-9, 1.4e-9, 1.4e-9, 1.4e-9)
#fm <-elast.comp.radius(tig.dtd.scale[1:4,"NEP-CCCTTAGGTCAAACGG"],S.radii,TRUE)
#
# compute radii for multiple cases
res<-elast.comp.radii(tig.dtd.scale[1:4,],S.radii,TRUE)
Res <- res[[1]]
dtd.predict<-res[[2]]
#
# integrate elast result into Seurat object
# 

source("elast.integrated.R")

#
# extract explanatory variables via glmnet
#
source("elast.glmnet.R")
en.model.nonzero.beta <- subset(en.model.beta, subset=abs(s0)>0.001)
#source("elast.biomaRt.R")
en.model.plus.beta <- subset(en.model.beta, subset=s0> 0.001)
en.model.minus.beta <- subset(en.model.beta, subset=s0< -0.001)

source(file.path(rdir,"elast.integrated.cor.visualize.R"))

source(file.path(rdir,"elast.integrated.cor.R"))

source(file.path(rdir,"elast.integrated.cor.visualize.R"))

source(file.path(rdir,"elast.integrated.cor.diffusionmap.R"))

 source("elast.integrate.data.R")
# tig.nep <- elast.integrate.data(tig.nep,Res)
# DimPlot(tig.nep,reduction = "pca")+FeaturePlot(tig.nep,features = "radii",
#                                                min.cutoff = 10,
#                                                max.cutoff = 60,reduction = "pca")
# VlnPlot(tig.nep,features = "radii")+scale_y_log10()+ylim(c(1,50))
# 
# VlnPlot(tig.nep,features = "pval")+ylim(c(0,0.2))#+scale_y_log10()
# RidgePlot(tig.nep,features="pval")
# FeaturePlot(tig.nep,features = "log.radii")

source("bulk.tig.mono.vs.co.R")
DimPlot(tig.bulk.seurat, reduction="pca", group.by = "celltype") + 
  FeaturePlot(tig.bulk.seurat,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"))+
  FeaturePlot(tig.nep,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"))


