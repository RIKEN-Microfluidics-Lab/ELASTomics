elast.integrate.data <- function(tig.nep,res){
  cellids <- matrix(nrow = ncol(tig.nep), ncol = 4)
  cellids <- data.frame(cellids)
  rownames(cellids)<-colnames(tig.nep)
  cellids[rownames(res) ,] <- res[,]
  colnames(cellids)<- c("rp","SE","t","pval")
  cellids[is.na(cellids)]<-0
  tig.nep[["radii"]]<-cellids$rp
  tig.nep[["se"]]<-cellids$SE
  tig.nep[["pval"]]<-cellids$pval
  return(tig.nep)
}
