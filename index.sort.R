datadir <- "/home/samba/sanger/shintaku/ELASTomics/20220922HiSeqX015_index"
index <- Read10X(data.dir = datadir)
index.dtd <- Read10X(data.dir = file.path(datadir,"CITE-seq"))

index.dtd.control<- index.dtd[,str_sub(colnames(index.dtd),1,10) == "Ind-DTDmix"]

index.dtd.subset <- index.dtd[,str_sub(colnames(index.dtd),9,21) %in% str_sub(colnames(index),9,21)]
index.subset <- index[,str_sub(colnames(index),9,21) %in% str_sub(colnames(index.dtd),9,21)]
colnames(index.subset) <- paste0("MF10Indc",str_sub(colnames(index.subset),9,21))
colnames(index.dtd.subset) <- paste0("MF10Indc",str_sub(colnames(index.dtd.subset),9,21))
index.seurat <- CreateSeuratObject(index.subset)
index.dtd.seurat <- CreateAssayObject(counts=index.dtd.subset)
index.seurat[["DTD"]]<-index.dtd.seurat

# annotate rt barcode 1-16
index.seurat[["barcode"]]<-str_sub(colnames(index.seurat),12,21)
rt.barcode <-data.frame(unlist(index.seurat[["barcode"]]))
colnames(rt.barcode)<-"represent"
source(file.path(rdir,"util/whitelist_encode.R"))
barcode<-unlist(read.table(file.path(rdir,"cell_id_list.txt")))
romin <- whitelist.umi_tools.encode(unique(rt.barcode$represent),barcode)
active_barcode <-barcode[sort(romin$index)]
index.seurat[["RT"]]<-whitelist.umi_tools.encode(rt.barcode$represent,active_barcode)$index
# annotate plate number A or B
plate.num <-data.frame(as.numeric(str_sub(colnames(index.seurat),9,10)))
colnames(plate.num)<-"number"
plate.num$plate <- "A"
plate.num[plate.num$number>6,]$plate <-"B"
plate.num$pool <- plate.num$number
plate.num[plate.num$number>6,]$pool <-plate.num[plate.num$number>6,]$number-6
index.seurat[["plate"]]<-plate.num$plate
index.seurat[["pool"]]<-plate.num$pool
cellids <-paste0(unlist(str_sub(colnames(index.seurat),1,8)),plate.num$plate,"-"
                 ,plate.num$pool,"-",unlist(index.seurat[["RT"]]))

index.seurat<-RenameCells(index.seurat,new.names = cellids)

source("/home/samba/public/shintaku/github/ELASTomics/io/hunter_Seurat_load_adt_data.R")
#load FACS index data file name is 9 char long.
indexdir <- "/home/samba/sanger/Shiomi/hunterindex"
files <- data.frame(list.files(file.path(datadir,"count"),pattern="counts.tsv.gz"))
colnames(files)<-"name"
batch <- unique(substr(files$name,1,8))
indexfiles <- file.path(indexdir,list.files(indexdir,pattern=batch))
channels <- rep(list(c("Events","FSC","SSC","Venus","PI","Hoechst")),length(indexfiles))

adt.csv <- cbind(load.adt(indexfiles[1],str_sub(list.files(indexdir,pattern=batch),1,9)[1],channels[[1]]),
                 load.adt(indexfiles[2],str_sub(list.files(indexdir,pattern=batch),1,9)[2],channels[[2]]))

adt.seurat <- CreateAssayObject(adt.csv[,colnames(adt.csv) %in% cellids])

index.seurat [["facs"]]<- adt.seurat
index.seurat <- NormalizeData(index.seurat,normalization.method ="CLR",assay = "DTD")
index.seurat <- NormalizeData(index.seurat,normalization.method ="CLR",assay = "facs")

# "FLD004-AACGTGAT" "FLD010-AAACATCG" "FLD040-ATGCCTAA" "FLD070-AGTGGTCA" "FLD150-ACCACTGT" "FLD500-ACATTGGC"
FeatureScatter(index.seurat,feature1="facs_Venus",feature2="FLD150-ACCACTGT")
VlnPlot(index.seurat,features=c("FLD004-AACGTGAT","FLD010-AAACATCG","FLD040-ATGCCTAA","FLD070-AGTGGTCA","FLD150-ACCACTGT","FLD500-ACATTGGC"),slot="counts")
