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
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 2000)
# normalize the dtd data with "RC" option.
tig <- NormalizeData(tig,assay="DTD",normalization.method = "RC",scale.factor = 1e2)
VlnPlot(tig,feature="dtd_FLD004",group.by = "condition")
#subset elastomics data after the normalization
tig.nep<-subset(tig,subset=condition=="NEP")
#
# re-scale the dtd data with the concentration of the dtd molecules in the solution
# 
source("elast.rescale.dtd.R")
concentration<- data.frame(c(2,6,9,3,9,3,1,1,1,1,1,1,1))
tig.nep<-rescale.dtd.data(tig.nep,concentration)
tig.dtd.scale <- tig.dtd.scale[c("FLD004","FLD010","FLD070","FLD500"),]
#
# compute a radius for a single case
source("elast.comp.radii.R")
# fm returns the nls output
# fm <-elast.comp.radius(tig.dtd.scale[1:4,1],TRUE)
#
# compute radii for multipe cases
# res returns a summary of the computed radii
res<-elast.comp.radii(tig.dtd.scale[1:4,],FALSE)
#
#cellids<-data.frame(colnames(tig.dtd.scale))
cellids <- matrix(nrow = ncol(tig.dtd.scale), ncol = 4)
cellids <- data.frame(cellids)
rownames(cellids)<-colnames(tig.dtd.scale)
cellids[rownames(res) ,] <- res[,]
colnames(cellids)<- c("rp","SE","t","pval")
cellids[is.na(cellids)]<-0

tig.nep[["radii"]]<-cellids$rp
tig.nep[["se"]]<-cellids$SE
tig.nep[["pval"]]<-cellids$pval


#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",margin=2,scale.factor = 1e5)
#pbmc <- ScaleData(pbmc,assay="ADT")

# annotate cells
source(file.path(rdir,"shiomi_Seurat_annotate_cells.R"))

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'shiomi_Seurat_10x_technical.R'))

# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/public/Shiomi/ELASTomicsindex"
channel <- c("Events","FSC","SSC","Venus","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")
source(file.path(rdir,'/io/hunter_Seurat_load_adt_data.R'))


#
# load cite-seq-count=FLD data as pbmc.tag
#
source(file.path(rdir,'io/hunter_Seurat_load_fld_data.R'))
# annotate cells with HTO and create FLD total
source(file.path(rdir,'20210816HiSeqX004_annotate_condition.R'))
#
#pbmc[["FLD"]] <- CreateAssayObject(counts=new_pbmc.tag[c("FLD004","FLD010","FLD040","FLD070","FLD150","Exter1","Exter2","Exter3","Exter4","Unmapp"),])
# check FLD by visualizing results
source(file.path(rdir,'20210816HiSeqX004_visualize_condition.R'))
# add HTO to Seurat object
pbmc[["HTO"]] <- CreateAssayObject(counts=pbmc.tag[c("T20CTL","T50CTL","TAZCTL"),])
#pbmc <- NormalizeData(pbmc, normalization.method = "CLR", scale.factor = 1e5,assay = "HTO")

DimPlot(pbmc,reduction="pca")
#p1<-DimPlot(pbmc,reduction="umap")
FeaturePlot(pbmc,features=c("FLD500","T20CTL","T50CTL","TAZCTL"),reduction = "pca")

#
source(file.path(rdir,"shiomi_fld_external_control_analysis.R"))
# first overview
# analyze data with PCA and UMAP
# find clusters and marker genes
#source(file.path(rdir,"shiomi_Seurat_technicalcheck.R"))

# check cell cycle dependence 
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
# find marker genes
source(file.path(rdir,"shiomi_Seurat_Marker.genes.R"))
# compute pseudotime and order cells along the gene expression
source(file.path(rdir,"shiomi_Seurat_monocle_pseudotime.R"))

#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))






