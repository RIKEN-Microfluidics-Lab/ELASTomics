library(stats)
library(sys)
# compute a single radius
elast.comp.radius <- function(DTD.scaled,S.radius,plot.flag){
  #DTD.scaled<-tig.dtd.scale[1:4,1]
  a1 <- -1.2167
  a2 <- 1.5336
  a3 <- -22.5083
  a4 <- -5.6117
  a5 <- -0.3363
  a6 <- -1.2160
  a7 <- 1.6470
  b1 <- 9 * sqrt(2) * pi * pi / 4
  fm<-nls(DTD.scaled ~ 6 * pi * DTD.scaled[1] * (1-(S.radius/P.radius))^2 /
           (b1 * (1-(S.radius/P.radius))^(-5/2) * (1 + a1 * (1-(S.radius/P.radius)) + a2 * (1-(S.radius/P.radius))^2) +
              a3 +
              a4 * (S.radius/P.radius) +
              a5 * (S.radius/P.radius)^2 +
              a6 * (S.radius/P.radius)^3 +
              a7 * (S.radius/P.radius)^4),start=c(P.radius=20),trace = FALSE)
  if (plot.flag){
    plot(S.radius,DTD.scaled,xlim=c(0,20),ylim=c(0,20))
    lines(S.radius,predict(fm))
    Sys.sleep(0.1)
  }
  return(fm)
}
#
# compute multipe radii
#
elast.comp.radii <- function(DTD.scaled,S.radii,plot.flag){
  #Calculate pore size
  cellids <- colnames(DTD.scaled)
  Res <- matrix(nrow = 0, ncol = 4)# 4 is for the four of outputs
  dtd.predict <- matrix(nrow=4,ncol=0)
  rownames(dtd.predict)<-rownames(DTD.scaled)
  errorid<-matrix(nrow=0,ncol=0)
  for (t in 1:ncol(DTD.scaled)){
    tryCatch({
      #Rat <- c(DTD.scaled[1:4,t])
        fm<-elast.comp.radius(c(DTD.scaled[,t]),S.radii,plot.flag)
        f.stat <- summary(fm)$coefficients
        fm.predict <- data.frame(unlist(predict(fm)))
        colnames(fm.predict)<-cellids[t]
        dtd.predict<-cbind(dtd.predict,fm.predict)
        rownames(f.stat) <- cellids[t]
        Res <- rbind(Res, f.stat)

    },
    error = function(e) {
      errorid<-c(errorid,cellids[t])
    }
    )
  }
  Res <- as.data.frame(Res)
  
  return(list(Res,dtd.predict))
}

