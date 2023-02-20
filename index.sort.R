library(stringdist)
rdir <- "/home/samba/public/shiomi/ELASTomics/"
datadir <- "/home/samba/sanger/shintaku/ELASTomics/20220922HiSeqX015_index/"
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

source(file.path(rdir,"io/hunter_Seurat_load_adt_data.R"))
#load FACS index data file name is 9 char long.
indexdir <- "/home/samba/sanger/Shiomi/hunterindex"
files <- data.frame(list.files(file.path(datadir,"count"),pattern="counts.tsv.gz"))
colnames(files)<-"name"
batch <- unique(substr(files$name,1,8))
indexfiles <- file.path(indexdir,list.files(indexdir,pattern=batch))
channels <- rep(list(c("Events","FSC","SSC","Venus","PI","Hoechst")),length(indexfiles))

adt.csv <- cbind(load.adt(indexfiles[1],str_sub(list.files(indexdir,pattern=batch),1,9)[1],channels[[1]]),
                 load.adt(indexfiles[2],str_sub(list.files(indexdir,pattern=batch),1,9)[2],channels[[2]]))
#adt.seurat <- CreateAssayObject(adt.csv[,colnames(adt.csv) %in% cellids])
adt.seurat <- as.matrix(adt.csv[,colnames(adt.csv) %in% cellids])
adt.seurat1 <- matrix(nrow = dim(adt.seurat)[1], ncol = dim(adt.seurat)[2])
for (icnt in 1:dim(adt.seurat)[1]){
  adt.seurat1[icnt,] <- as.integer(adt.seurat[icnt,])
}
colnames(adt.seurat1) <- colnames(adt.seurat)
rownames(adt.seurat1) <- rownames(adt.seurat)
adt.seurat <- CreateAssayObject(adt.seurat1)

index.seurat [["facs"]]<- adt.seurat
index.seurat <- NormalizeData(index.seurat,normalization.method ="CLR",assay = "DTD")
index.seurat <- NormalizeData(index.seurat,normalization.method ="CLR",assay = "facs")

# "FLD004-AACGTGAT" "FLD010-AAACATCG" "FLD040-ATGCCTAA" "FLD070-AGTGGTCA" "FLD150-ACCACTGT" "FLD500-ACATTGGC"
FeatureScatter(index.seurat,feature1="facs_Venus",feature2="FLD500-ACATTGGC")
VlnPlot(index.seurat,features=c("FLD004-AACGTGAT","FLD010-AAACATCG","FLD040-ATGCCTAA","FLD070-AGTGGTCA","FLD150-ACCACTGT","FLD500-ACATTGGC"),slot="counts")

#DTDtable <- bind_cols(as.data.frame(t(index.seurat[["DTD"]][c("FLD004-AACGTGAT", "FLD010-AAACATCG", "FLD040-ATGCCTAA", "FLD070-AGTGGTCA", "FLD150-ACCACTGT", "FLD500-ACATTGGC")])), as.data.frame(t(index.seurat[["facs"]]["Venus"])))
DTDtable <- bind_cols(as.data.frame(t(index.seurat[["DTD"]]@counts)), as.data.frame(t(index.seurat[["facs"]]@counts)))
DTDtable <- DTDtable[c("FLD004-AACGTGAT", "FLD010-AAACATCG", "FLD040-ATGCCTAA", "FLD070-AGTGGTCA", "FLD150-ACCACTGT", "FLD500-ACATTGGC", "Venus")]
colnames(DTDtable) <- c("FLD004","FLD010","FLD040","FLD070","FLD150","FLD500","Venus")
DTDtable$Sum <- DTDtable$FLD004+DTDtable$FLD010+DTDtable$FLD040+DTDtable$FLD070+DTDtable$FLD150+DTDtable$FLD500
DTDtable$Plot <- str_sub(rownames(DTDtable), start = 9, end = 9)
DTDtable$Plot1 <- if_else(as.integer(str_sub(rownames(DTDtable), start = 11, end = 11)) < 4, true = 1, false = 2)
DTDtable$Plot <- paste(DTDtable$Plot, DTDtable$Plot1, sep = "") 

g <- ggplot(DTDtable, aes(y = FLD500, x = Plot),color = Plot) + geom_point(position = "identity", alpha = 0.8) +
  geom_boxplot(color = c("#EA5514", "#F39800", "#00A0E9", "#036EB8"), alpha = 1) +
  geom_jitter(size = 2, width = 0.2) + scale_y_log10()+theme_classic()
plot(g)

g <- ggplot(DTDtable, aes(y = FLD500, x = Venus, color = Plot))+geom_point()+ scale_color_manual(values = c("#EA5514", "#F39800", "#00A0E9", "#036EB8"))+theme_bw()
plot(g)  

DTDtable %>%
  ggplot(aes(x = Venus, y = FLD500)) +
  geom_point() +
  geom_density_2d()+scale_y_log10()+scale_x_log10()+theme_classic()
cor(DTDtable$Venus, DTDtable$FLD500, method="spearman")

+scale_y_sqrt()+theme_bw()
  theme_minimal(base_size = 20) 

  