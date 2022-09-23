load.elast.data <- function(wdir,cellname,min.cell.num){
  data <- Read10X(data.dir = wdir)
  tig <- CreateSeuratObject(data$`Gene Expression`,min.cells = min.cell.num)
  if (!is_empty(data$Custom)){
    tig[["DTD"]]<- CreateAssayObject(counts = data$Custom)
  }
  if (!is_empty(data$`Antibody Capture`)){
    tig[["ADT"]]<- CreateAssayObject(counts = data$`Antibody Capture`)
  }
  cellids <-colnames(tig)
  cellids<-paste0(cellname,substr(cellids,1,16))
  tig<-RenameCells(tig, new.names = cellids)
  return(tig)
}
