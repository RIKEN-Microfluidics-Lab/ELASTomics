rescale.dtd.data <- function(tig.nep,concentration){
  tig.dtd.data <- data.frame(tig.nep[["DTD"]]@data)
  rownames(concentration)<-rownames(tig.dtd.data)
  colnames(concentration)<-"concentration"
  tig.dtd.scale<-sweep(as.matrix(tig.dtd.data),1,as.matrix(concentration),"/")
  colnames(tig.dtd.scale)<-colnames(tig.nep)
  tig.nep[["DTD"]]@data<-tig.dtd.scale
  return(tig.nep)
}
