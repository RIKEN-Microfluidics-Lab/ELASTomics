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
nep <-load.elast.data(wdir,"NEP-",100)
#nep[["condition"]]<-"NEP"
#
nep <- NormalizeData(nep, normalization.method = "LogNormalize", scale.factor = 1e5)
# mitochondrial gene percent and remove dead cells
nep[["percent.mt"]] <- PercentageFeatureSet(nep, pattern = "^mt-")
nep <- subset(nep, subset= percent.mt<5)
nep <- subset(nep,subset=nCount_RNA >1000)
#
#
# merge all the tig data
#
nep <- FindVariableFeatures(nep, selection.method = "vst", nfeatures = 300)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nep), 40)
# plot variable features with labels
plot1 <- VariableFeaturePlot(nep)
plot1 <- LabelPoints(plot = nep, points = top10)
plot1
# normalize the dtd data with "RC" option.
nep <- NormalizeData(nep,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(nep)
nep <- ScaleData(nep, features = all.genes)
nep <- RunPCA(nep, npcs=20, features = VariableFeatures(object = nep))
DimPlot(nep, reduction = "pca",group.by = "condition")
#source("Seurat.clustering.R")
JackStrawPlot(nep, dims = 1:20)
ElbowPlot(nep)

nep <- FindNeighbors(nep, dims = 1:18)
nep <- FindClusters(nep, resolution = 0.3)

nep <- RunUMAP(nep, dims = 1:19)
FeaturePlot(nep,features = c("Plek","Cebpa","Epor","Hlf","Mpo","Spib","Cd79a","Dntt","Mpl","Pf4","Vwf","Itga2b","dtd_ADT130","dtd_FLD004","nFeature_RNA"),reduction="umap")+DimPlot(nep)
FeaturePlot(nep,features="dtd_FLD004",max.cutoff = 2)

#subset elastomics data after the normalization
tig.nep<-subset(tig,subset=condition=="NEP")
tig.ctl<-subset(tig,subset=condition=="CTL")

#
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

# source("elast.integrate.data.R")
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


