library(ggplot2)
library(stringr)
library(reshape2)
library(gridExtra)
# create a list of whitelist files
filelist_whitelist <- data.frame(list.files(datadir,pattern="whitelist.txt"))
colnames(filelist_whitelist) <- c("filename")
# read whitelist file
for (icnt in 1:nrow(filelist_whitelist)){
  whitelist <- read.table(file.path(datadir,filelist_whitelist[icnt,]), sep="\t")
  colnames(whitelist) <- c('represent','variants', 'total','counts')
  # create functions for encoding
  encoded <- whitelist.umi_tools.list(whitelist,barcode$V1)
  encoded <- encoded[order(encoded$first_index,decreasing = FALSE),]
  correct_encoded <- encoded[encoded$lv_total<1,]
  correct_encoded$batch <- str_sub(filelist_whitelist[icnt,],1,10)
  rownames(correct_encoded) <- paste(correct_encoded$batch,correct_encoded$first_index,sep="-")
  if (icnt>1){
    allencoded <- rbind(allencoded,correct_encoded)
  }
  else{
    allencoded <- correct_encoded
  }
}


allencoded$GC <- as.numeric(lapply(lapply(as.character(allencoded$first_barcode),s2c),GC))
p0 <- ggplot(allencoded,aes(x=factor(GC),y=count))+ geom_violin()+ scale_y_log10()+geom_boxplot()
p1 <- ggplot(allencoded, aes(x = first_index, y = count, color=first_barcode))+geom_violin()
p2 <- ggplot(allencoded, aes(x = batch, y = count, color=batch))+geom_violin()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

grid.arrange(p0,p1, p2, nrow = 3)
rm(correct_encoded,encoded,p0,p1,p2)
