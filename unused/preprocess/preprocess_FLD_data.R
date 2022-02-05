load.fld <- function(datadir,count_type,barcode,decode){
  sampleID <- list.dirs(path=file.path(datadir,"CITE-seq"), full.names = FALSE, recursive = FALSE)
  for (i in 1:length(sampleID)){
    #FLDmap <- matrix(nrow=7, ncol=0)
    #count_type <- "read_count"
    samfol = file.path(datadir, "CITE-seq", sampleID[i], count_type,sep="")
    FLD.data <- Read10X(data.dir = samfol, gene.column=1)
    FLData <- as.data.frame(as.matrix(FLD.data))
    if (decode==TRUE){
      romin <- whitelist.umi_tools.encode(colnames(FLD.data),barcode$V1)
      romin1 <- romin$index * ifelse(romin$value>2, 0, 1)
      FLData <- as.data.frame(cbind(t(FLData), romin1))
      FLData = FLData[FLData$romin1!=0,]
    
      FLData <- group_by(FLData, romin1)
      FLData <- summarise_each(FLData, dplyr::funs(sum))
      FLDmap <- as.data.frame(FLData)
      FLDmap$GC <- barcode[FLDmap$romin1,]$GC
      iName = substr(sampleID[i], 1, 10)
      rownames(FLDmap) <- paste0(iName,"-",FLDmap$romin1 )
    }else{
      FLDmap <- as.data.frame(FLData)
    }
    
    
    if (i==1){
      FLDmapALL<-FLDmap
    }else{
      FLDmapALL <-rbind(FLDmapALL, FLDmap)
    }
  }
  
  FLDmapALL <- as.data.frame(FLDmapALL)
  return(FLDmapALL)
}

