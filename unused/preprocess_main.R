
# decode the single cell data from whitelist of UMI-tools output
#datadir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_TIG/"
#datadir <- "/home/samba/sanger/Shiomi/20211124HiSeqX006_TIG/"

#datadir <-"/home/samba/sanger/shintaku/20220109HiSeqX008_tig/"
#wdir <- "/home/samba/sanger/shintaku/20210323MiSeq015Ana10X/"
#wdir <- "/home/samba/public/shintaku/ELASTomics/"
#rdir <- "/home/samba/public/shintaku/github/ELASTomics/"

#
# cell_id_list2.txt contains all barcodes
# cell_id_list.txt contains selected barcodes by GC percent.
#barcode <- read.table(file.path(rdir,"cell_id_list2.txt"))
#barcode <- read.table(file.path("/home/samba/sanger/Program/cellranger-6.1.2/lib/python/cellranger/barcodes/3M-february-2018.txt"))
#rownames(barcode)<-barcode$V1
#barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

# 10x for big data, facs for small data
#source(file.path(rdir,"10x_first_data_process.R"))
# preprocess FLD data
#source(file.path(rdir,"preprocess/preprocess_FLD_data.R"))

#HiSeq006

#sampleID <- list.dirs(path=file.path(datadir,"CITE-seq"), full.names = FALSE, recursive = FALSE)
# #count_type <- "read_count"
#count_type <- "umi_count"
#FLDmapALL <- load.fld(datadir,count_type,barcode,FALSE)
# 
# samfol = file.path(datadir, "CITE-seq", sampleID[4], count_type,sep="")
# FLD.data <- Read10X(data.dir = samfol, gene.column=1)
# FLData<- as.data.frame(as.matrix(FLD.data))
# FLData<- FLData[1:6, ]
#samfol = file.path(datadir2, "CITE-seq", "AS-ExtC-2E", count_type,sep="")
#FLD.data <- Read10X(data.dir = samfol, gene.column=1)
#FLData1<- as.data.frame(as.matrix(FLD.data))
#FLData1<- FLData1[9:12, ]
#FLData <- bind_rows(FLData, FLData1)
#rm(FLData1)
#FLDmap <- as.data.frame(FLData)
#colnames(FLDmap) <- paste0("C02-TIG-EP","-",colnames(FLDmap))

#samfol = file.path(datadir1, "CITE-seq", sampleID[2], count_type,sep="")
#FLD.data <- Read10X(data.dir = samfol, gene.column=1)
#FLData<- as.data.frame(as.matrix(FLD.data))
#FLData<- FLData[1:6, ]
#samfol = file.path(datadir2, "CITE-seq", "AS-ExtC-1C", count_type,sep="")
#FLD.data <- Read10X(data.dir = samfol, gene.column=1)
#FLData1<- as.data.frame(as.matrix(FLD.data))
#FLData1<- FLData1[9:12, ]
#FLData <- bind_rows(FLData, FLData1)
#rm(FLData1)
#FLDmap0 <- as.data.frame(FLData)
#colnames(FLDmap0) <- paste0("C02-TIG-C.","-",colnames(FLDmap0))
#FLDmapALL <- bind_cols(FLDmap0, FLDmap)
#rm(FLData, FLDmap, FLDmap0)

pbmctag <- matrix(nrow = 10, ncol = length(cellids))
rownames(pbmctag) <- c("FLD004", "FLD010", "FLD040", "FLD070", "FLD150", "FLD500", "EX1CTL",  "EX2CTL", "EX3CTL",  "EX4CTL")
colnames(pbmctag) <- cellids
for (i in 1:length(cellids)){
  tryCatch(
    {pbmctag[1:10,i] <- FLDmapALL[, cellids[i]]}
    , error = function(e) {print(i)}
  )
}
rm(FLDmapALL)

tig[["DTD"]]

# what are these numbers? concentrations
pbmctagR <-as.data.frame(t(pbmctag))
pbmctagR$FLD004 <- (1000 * pbmctagR$FLD004 / 200)
pbmctagR$FLD010 <- (1000 * pbmctagR$FLD010 / 600)
pbmctagR$FLD040 <- (1000 * pbmctagR$FLD070 / 900)
pbmctagR$FLD070 <- (1000 * pbmctagR$FLD500 / 300)
pbmctagR$FLD150 <- (1000 * pbmctagR$FLD150 / 900)
pbmctagR$FLD500 <- (1000 * pbmctagR$FLD500 / 300)

pbmctagR$FLDALL <- apply(pbmctagR[1:6], 1, sum)
pbmctagR$pFLD004 <- log10((pbmctagR$FLD004)+1)
pbmctagR$pFLD500 <- log10((pbmctagR$FLD500)+1)
pbmctagR$pFLDRatio <- 100 *(pbmctagR$FLD500)/((pbmctagR$FLD004)+1)
pbmctagR$pFLDALL <- log10((pbmctagR$FLDALL)+1)

#normalize by total DTD
pbmctagR$RFLD004 <- (100 * pbmctagR$FLD004 / ((pbmctagR$FLDALL)+1))
pbmctagR$RFLD010 <- (100 * pbmctagR$FLD010 / ((pbmctagR$FLDALL)+1))
pbmctagR$RFLD040 <- (100 * pbmctagR$FLD040 / ((pbmctagR$FLDALL)+1))
pbmctagR$RFLD070 <- (100 * pbmctagR$FLD070 / ((pbmctagR$FLDALL)+1))
pbmctagR$RFLD150 <- (100 * pbmctagR$FLD150 / ((pbmctagR$FLDALL)+1))
pbmctagR$RFLD500 <- (100 * pbmctagR$FLD500 / ((pbmctagR$FLDALL)+1))

# take the elp part
pbmctagE <-  pbmctagR[5292:10383, ]
pbmctagE <-  filter(pbmctagE, dplyr::between(pFLDRatio, 0, 200))
pbmctagE <- pbmctagE[order(pbmctagE$pFLDRatio),]
pbmctagMean <- matrix(nrow = 0, ncol = 22)
for (p in 1:10) {
  pbmctagMean <- rbind(pbmctagMean, colMeans(pbmctagE[(500 * (p-1)):(500 * p),]))
}
pbmctagMean <- as.data.frame(pbmctagMean)
pbmctagMean <- pbmctagMean[, 14:21]
colnames(pbmctagMean) <-colnames(pbmctagE)[14:21]
pbmctagMean <- pbmctagMean[, c(3,4,6,8)]

g <- ggplot(pbmctagE, aes(x = pFLDALL, y = pFLDRatio))
g <- g + geom_point() +ylim(0,60)
plot(g)

#
pbmctagE[, 16:21] %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=Val, fill=MesureType)) +
  geom_violin() #+ scale_y_log10() 


#Calculate pore size
Dex <- c(1.4, 2.7, 6.3, 15.1)
a1 <- -1.2167
a2 <- 1.5336
a3 <- -22.5083
a4 <- -5.6117
a5 <- -0.3363
a6 <- -1.2160
a7 <- 1.6470
b1 <- 9 * sqrt(2) * pi * pi / 4
Res <- matrix(nrow = 0, ncol = 4)
for (t in 1:nrow(pbmctagE)){
  tryCatch({
    Rat <- c(pbmctagE[t, 16], pbmctagE[t, 17], pbmctagE[t, 19], pbmctagE[t, 21])
    fm<-nls(Rat ~ 6 * pi * Rat[1] * (1-(Dex/Por))^2 / (b1 * (1-(Dex/Por))^(-5/2) * (1 + a1 * (1-(Dex/Por)) + a2 * (1-(Dex/Por))^2) + a3 + a4 * (Dex/Por) + a5 * (Dex/Por)^2 + a6 * (Dex/Por)^3 + a7 * (Dex/Por)^4),start=c(Por=1000),trace = FALSE)
    if (0 < coef(fm) && 100 > coef(fm)) {
      f.stat <- summary(fm)$coefficients
      rownames(f.stat) <- rownames(pbmctagE)[t]
      Res <- rbind(Res, f.stat)
    } 
  }
  , error = function(e) {print(t)}
  )
}
Res <- as.data.frame(Res)
colnames(Res) <-c("Estimate", "Std", "tvalue", "pvalue")
g <- ggplot(Res, aes(x = pvalue))
g <- g + geom_histogram()#+xlim(0,60)
plot(g)

g <- ggplot(Res, aes(x = Estimate, y = pvalue))
g <- g + geom_point() #+ylim(0,60)
plot(g)



pbmctagR <-as.data.frame(t(pbmctagR))
pbmc[["FLD"]] <- CreateAssayObject(counts=pbmctagR)
pbmc[['type']] <- substr(cellids,9,9)

pbmctagR <-as.data.frame(t(pbmctagR))
pbmctagR$Amu <- ifelse(pbmctagR$pFLDALL<0.1,  "None",ifelse(pbmctagR$pFLDALL<4, "EP", "Hole"))

pbmc[['Amu']] <- pbmctagR$Amu
