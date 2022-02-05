library(stats)
library(sys)
elast.comp.radius <- function(DTD.scaled,plot.flag){
  S.radius <- c(1.4, 2.7, 6.3, 15.1) # Stokes radii of DTD
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
              a7 * (S.radius/P.radius)^4),start=c(P.radius=1000),trace = FALSE)
  if (plot.flag){
    plot(S.radius,DTD.scaled,xlim=c(0,15),ylim=c(0,20))
    lines(S.radius,predict(fm))
    Sys.sleep(0.1)
  }
  return(fm)
}
elast.comp.radii <- function(DTD.scaled,plot.flag){
  #Calculate pore size
  cellids <- colnames(DTD.scaled)
  Res <- matrix(nrow = 0, ncol = 4)
  errorid<-matrix(nrow=0,ncol=0)
  for (t in 1:ncol(DTD.scaled)){
    tryCatch({
      #Rat <- c(DTD.scaled[1:4,t])
        fm<-elast.comp.radius(c(DTD.scaled[1:4,t]),plot.flag)
        f.stat <- summary(fm)$coefficients
        rownames(f.stat) <- cellids[t]
        Res <- rbind(Res, f.stat)

    },
    error = function(e) {
      errorid<-c(errorid,cellids[t])
    }
    )
  }
  Res <- as.data.frame(Res)

  return(Res)

}

