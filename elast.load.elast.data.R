load.elast.data <- function(wdir,cellname){
  data <- Read10X(data.dir = wdir)
  tig <- CreateSeuratObject(data$`Gene Expression`)
  tig[["DTD"]]<- CreateAssayObject(counts = data$Custom)
  cellids <-colnames(tig)
  cellids<-paste0(cellname,substr(cellids,1,16))
  tig<-RenameCells(tig, new.names = cellids)
  return(tig)
}