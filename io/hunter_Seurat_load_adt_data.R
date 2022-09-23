load.adt <- function (indexfiles,batch,channel) {
  rowend <- length(channel)+1
  adt.xlsx <- read.xlsx(indexfiles[1],sheet=1, rows=16:111,cols=1:rowend,colNames = FALSE,rowNames = TRUE)
  colnames(adt.xlsx)<- channel#c("Events","FSC","SSC","Venus","APC","mCherry")
  
  # create normalized GFP
  #adt.xlsx$normGFP <- (adt.xlsx$Venus/adt.xlsx$mCherry)
  
  
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


