require(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils)
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
wdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_TIG/CNTRL/outs/filtered_feature_bc_matrix/"
ctl <-load.elast.data(wdir,"CTL-")
ctl[["condition"]]<-"CTL"
# elastomics data
wdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_TIG/EP/outs/filtered_feature_bc_matrix/"
nep <-load.elast.data(wdir,"NEP-")
nep[["condition"]]<-"NEP"
#
# merge all the tig data
#
tig <- merge(ctl, y=nep)
#
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)
# mitochondrial gene percent and remove dead cells
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^MT-")
tig <- subset(tig, subset= percent.mt<5)
#
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 5000)
# normalize the dtd data with "RC" option.
tig <- NormalizeData(tig,assay="DTD",normalization.method = "RC",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca",group.by = "condition")

#subset elastomics data after the normalization
tig.nep<-subset(tig,subset=condition=="NEP")
#
# re-scale the dtd data with the concentration of the dtd molecules in the solution
# 
source("elast.rescale.dtd.R")
# FLD004,FLD010,FLD040,FLD070,FLD150,FLD500, and others
concentration<- data.frame(c(2,6,9,3,9,3,1,1,1,1,1,1,1))
tig.nep<-rescale.dtd.data(tig.nep,concentration)
tig.dtdd.scale<-tig.nep[["DTD"]]@data
tig.dtd.scale <- tig.dtd.scale[c("FLD004","FLD010","FLD070","FLD500"),]
#
# compute a radius for a single case
source("elast.comp.radii.R")
#
# S.radii: input Stokes radii of DTD
# tig.dtd.scale; scaled amount of DTD imported to cells
# The nrow of the tig.dtd.scale must match with the length of S.radii
#
S.radii <- c(1.4, 2.7, 6.3, 15.1) # Stokes radii of DTD
# fm <-elast.comp.radius(tig.dtd.scale[1:4,1],S.radii,TRUE)
#
# compute radii for multipe cases
res<-elast.comp.radii(tig.dtd.scale[1:4,],S.radii,FALSE)
#
# integrate elast result into Seurat object
# 
tig.nep <- elast.integrate.data(tig.nep,res)
VlnPlot(tig.nep,features = "pval")+ylim(c(0,0.2))#+scale_y_log10()
RidgePlot(tig.nep,features="pval")
FeaturePlot(tig.nep,features = "log.radii")

source("elast.bulk.tig.mono.vs.co.R")
DimPlot(tig.bulk.seurat, reduction="pca", group.by = "celltype") + 
  FeaturePlot(tig.bulk.seurat,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"))+
  FeaturePlot(tig.nep,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","RAD51AP1"))




