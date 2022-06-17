require(readr)
library(plyr,dplyr,tidyr)
library(ggplot2)
library(tidyverse,R.utils, RCurl)
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
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/CNTRL/outs/filtered_feature_bc_matrix/"
ctl <-load.elast.data(wdir,"CTL-",100)
ctl[["condition"]]<-"CTL"
# elastomics data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20211124HiSeqX006_TIG/EP/outs/filtered_feature_bc_matrix/"
nep <-load.elast.data(wdir,"NEP-",100)
nep[["condition"]]<-"NEP"
#
tig <- merge(ctl, y=nep)
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)
# mitochondrial gene percent and remove dead cells
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^MT-")
tig <- subset(tig, subset= percent.mt<5)
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
# normalize the dtd data with "RC" option.
tig <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca",group.by = "condition")
source("Seurat.clustering.R")

#subset elastomics data after the normalization
tig.nep<-subset(tig,subset=condition=="NEP")
tig.ctl<-subset(tig,subset=condition=="CTL")
DimPlot(tig.nep, label = TRUE)
#
# extract explanatory variables via glmnet
#
source("elast.glmnet.R")
en.model.nonzero.beta <- subset(en.model.beta, subset=abs(s0)>0.001)
#source("elast.biomaRt.R")
en.model.plus.beta <- subset(en.model.beta, subset=s0> 0.001)
en.model.minus.beta <- subset(en.model.beta, subset=s0< -0.001)


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
S.radii <- c(1.4, 2.7, 6.3, 15.1) # Stokes radii of DTD
#fm <-elast.comp.radius(tig.dtd.scale[1:4,"NEP-CCCTTAGGTCAAACGG"],S.radii,TRUE)
#
# compute radii for multiple cases
res<-elast.comp.radii(tig.dtd.scale[1:4,],S.radii,TRUE)
Res <- res[[1]]
dtd.predict<-res[[2]]
#
# integrate elast result into Seurat object
# 

source("elast.integrated.R.R")

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


