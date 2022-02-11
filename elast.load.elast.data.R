load.elast.data <- function(wdir,cellname,min.cell.num){
  data <- Read10X(data.dir = wdir)
  tig <- CreateSeuratObject(data$`Gene Expression`,min.cells = min.cell.num)
  tig[["DTD"]]<- CreateAssayObject(counts = data$Custom)
  cellids <-colnames(tig)
  cellids<-paste0(cellname,substr(cellids,1,16))
  tig<-RenameCells(tig, new.names = cellids)
  return(tig)
}
