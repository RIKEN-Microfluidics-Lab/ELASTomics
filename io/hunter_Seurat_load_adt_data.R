load.adt <- function (indexfiles,batch,channel) {
  rowend <- length(channel)+1
  adt.xlsx <- read.xlsx(indexfiles[1],sheet=1, rows=16:111,cols=1:rowend,colNames = FALSE,rowNames = TRUE)
  colnames(adt.xlsx)<- channel#c("Events","FSC","SSC","Venus","APC","mCherry")
  
  # create normalized GFP
  adt.xlsx$normGFP <- (adt.xlsx$Venus/adt.xlsx$mCherry)
  
  
  wellid <- rownames(adt.xlsx)
  for (ij in 1:length(wellid)){
    icnt=as.numeric(which(LETTERS==substr(wellid[ij],1,1)))
    jcnt=as.numeric(substr(wellid[ij],2,3))
    pindex= (jcnt+1)%/%2
    rtindex <- icnt+((jcnt-1)%%2)*8
    adt.xlsx$cell[ij] <- paste0(batch,"-",pindex,"-",rtindex)
  }
  rownames(adt.xlsx)<-adt.xlsx$cell
  adt.xlsx <- t(adt.xlsx)
  adt.xlsx <- adt.xlsx[-nrow(adt.xlsx),]
  return(adt.xlsx)
}

#load FACS index data file name is 8 char long.
files <- data.frame(list.files(file.path(datadir,"count"),pattern="counts.tsv.gz"))
colnames(files)<-"name"
batch <- unique(substr(files$name,1,8))
indexfiles <- file.path(indexdir,paste0(batch,".xlsx"))
channels <- rep(list(channel),length(indexfiles))
adt.xlsx <- mapply(load.adt,indexfiles,batch,channels)
adt.csv <- do.call("cbind",adt.xlsx)

# adding adt data
seladt.csv <- adt.csv[,cellids] # extract cells detected in RNA-seq
seladt.csv[is.na(seladt.csv)] <- 0
pbmc.adt <-as.sparse(seladt.csv) # convert the format to sparse matrix
pbmc[["ADT"]] <- CreateAssayObject(counts=pbmc.adt)

#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",scale.factor = 1e5)
#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",margin=2,scale.factor = 1e5)
#NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)

#pbmc <- ScaleData(pbmc,assay="ADT")

rm(batch,indexfiles,channels,adt.xlsx,adt.csv,seladt.csv,pbmc.adt,channel)
