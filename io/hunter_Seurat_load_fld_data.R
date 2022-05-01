library(Seurat)
library(stringr)
library(dplyr)

#source(file.path(rdir,"preprocess/preprocess_FLD_data.R"))
decode=TRUE
FLDmapALL <- load.fld(datadir,"umi_count",active_barcode,decode)
sampleID <- list.dirs(path=file.path(datadir,"CITE-seq"), full.names = FALSE, recursive = FALSE)
samfol = file.path(datadir, "CITE-seq", sampleID[1], count_type,sep="")
FLD.data <- Read10X(data.dir = samfol, gene.column=1)
DTD <- CreateSeuratObject(FLD.data,assay="DTD")


tagnames <- colnames(FLDmapALL)
tagnames <- data.frame(strsplit(tagnames,"-"))
colnames(FLDmapALL) <- tagnames[1,]
rownames(FLDmapALL) <- toupper(rownames(FLDmapALL))
FLDmapALL<-t(FLDmapALL)
# extract cellids shared with cDNA
cellids_fld <- colnames(FLDmapALL)
selcellids <- intersect(cellids,cellids_fld)
seladt.FLD.csv <- FLDmapALL[,selcellids] # extract cells detected in RNA-seq
rm(FLDmapALL)
#create empty data frame with matchig rows and cols of RNA-seq
pbmc.tag<-data.frame((matrix(0,nrow=nrow(seladt.FLD.csv),ncol=length(cellids))))
rownames(pbmc.tag)<-rownames(seladt.FLD.csv)
colnames(pbmc.tag)<-cellids
#replace the created data frame with FLD data
pbmc.tag[,selcellids]<-seladt.FLD.csv[,selcellids]
rm(cellids_fld,selcellids)
#seladt.FLD.csv[is.na(seladt.FLD.csv$romin),]<-0



rm(seladt.FLD.csv,pbmc.adt.FLD,FLDcomb,FLDmelt)

