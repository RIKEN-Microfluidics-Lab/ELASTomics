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
wdir <-"/home/samba/sanger/Shiomi/20210601_HiSeqAna/"
datadir<-"/home/samba/sanger/Shiomi/20210601_HiSeq/"
barcode <- read.table(file.path("/home/samba/sanger/shintaku/github/hunter/cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'unused/preprocess/preprocess_whitelist.R'))
active_barcode <- barcode[sort(unique(allencoded$first_index)),]
encode_barcode=TRUE
index.sort <- Read10X(data.dir = wdir)
index.sort <-CreateSeuratObject(index.sort)
cellids <- colnames(index.sort)
index.sort[['plate']] <- substr(cellids,1,3)
index.sort[['condition']] <- substr(cellids,5,5)
index.sort[['cell']] <- substr(cellids,4,6)
index.sort[['gate']] <- substr(cellids,7,8)
index.sort[['pool']] <- substr(cellids,10,10)
index.sort[['rtid']] <- substr(cellids,12,13)

indexdir <- "/home/samba/sanger/Shiomi/hunterindex"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"io/hunter_Seurat_load_adt_data.R"))
index.sort[["ADT"]] <- CreateAssayObject(counts=pbmc.adt)

source(file.path(rdir,"unused/preprocess/preprocess_FLD_data.R"))
count_type="read_count"
source(file.path(rdir,"io/hunter_Seurat_load_fld_data.R"))
index.sort[["DTD"]]<-CreateAssayObject(counts=pbmc.tag)
index.sort<-NormalizeData(index.sort,assay="DTD",normalization.method = "LogNormalize",scale.factor = 1e2)
index.sort<-NormalizeData(index.sort,assay="ADT",normalization.method = "LogNormalize",scale.factor = 1e2)
FeatureScatter(subset(index.sort,subset=condition=="A",invert=TRUE),
               feature2="FLD004",feature1 = "Venus",group.by = "cell")+
FeatureScatter(index.sort,
               feature2="FLD004",feature1 = "Venus",group.by = "cell")
