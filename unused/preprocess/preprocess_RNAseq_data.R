
# load count data
#datadir <- dir #"/home/watson/sanger/shintaku/20210216HiSeqX002/count/"
files <- data.frame(list.files(file.path(datadir,"count"),pattern="counts.tsv.gz"))
colnames(files)<-"name"

for (icnt in 1:nrow(files)){
  #icnt=1
  myData <- read.table(file.path(datadir,"count",files[icnt,]), header = TRUE)
  if (encode_barcode==TRUE){
    white<-subset(allencoded,allencoded$batch==str_sub(files[icnt,],1,10))
    encoded <- whitelist.umi_tools.encode(myData$cell,active_barcode$V1)
    myData$cell <- paste0(str_sub(files[icnt,],1,10),"-",encoded$index)
    myData <- myData[encoded$value<1,]
  }else{
    myData$cell <- paste0(str_sub(files[icnt,],1,10),"-",myData$cell)
  }
  
  if (icnt>1){
    allData <- rbind(allData,myData)
  }
  else{
    allData <- myData
  }
}
rm(myData,white,encoded)
