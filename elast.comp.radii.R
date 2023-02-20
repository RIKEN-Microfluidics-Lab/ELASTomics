library(stats)
library(sys)
library(matlib)
library(numbers)
library(reshape2)
# compute a single radius
elast.comp.radius <- function(DTD.scaled,S.radii,par,alpha,scale,plot.flag){
  #DTD.scaled<-tig.dtd.scale[1:4,1]
  a1 <- -1.2167
  a2 <- 1.5336
  a3 <- -22.5083
  a4 <- -5.6117
  a5 <- -0.3363
  a6 <- -1.2160
  a7 <- 1.6470
  b1 <- 9 * sqrt(2) * pi * pi / 4
  C <-1.39e-46
  kbT<-300*1.380649e-23
  #S.radius<-S.radii
  #P.radius<-10e-9
  gamma<-par$gamma #1.138212e-13
  sigma<-par$sigma #5.359301e-06
  
  iter_max <-1000000
  delta_history<-matrix(nrow=101,ncol=3)
  par_history<-matrix(nrow=101,ncol=4)
  iter<-0
  eps <-1
  while ((iter < iter_max) & (eps > 1e-2)){
    Imatrix<-matrix(nrow = 5,ncol = length(S.radii))
    for (j in 1:length(S.radii)) {
      S.radius<-S.radii[j]
      Hp2<-function(P.radius) {P.radius^2*(6 * pi * (1-(S.radius/P.radius))^2 /
                                             (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                                (1 + a1 * (1-(S.radius/P.radius)) +
                                                   a2 * (1-(S.radius/P.radius))^2) +
                                                a3 +
                                                a4 * (S.radius/P.radius) +
                                                a5 * (S.radius/P.radius)^2 +
                                                a6 * (S.radius/P.radius)^3 +
                                                a7 * (S.radius/P.radius)^4))*
          exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-sigma*pi*P.radius^2)/kbT)}
      
      Hp3<-function(P.radius) {P.radius^3*(6 * pi * (1-(S.radius/P.radius))^2 /
                                             (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                                (1 + a1 * (1-(S.radius/P.radius)) +
                                                   a2 * (1-(S.radius/P.radius))^2) +
                                                a3 +
                                                a4 * (S.radius/P.radius) +
                                                a5 * (S.radius/P.radius)^2 +
                                                a6 * (S.radius/P.radius)^3 +
                                                a7 * (S.radius/P.radius)^4))*
          exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-sigma*pi*P.radius^2)/kbT)*2*pi/kbT}

      Hp4<-function(P.radius) {P.radius^4*(6 * pi * (1-(S.radius/P.radius))^2 /
                                             (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                                (1 + a1 * (1-(S.radius/P.radius)) + 
                                                   a2 * (1-(S.radius/P.radius))^2) +
                                                a3 +
                                                a4 * (S.radius/P.radius) +
                                                a5 * (S.radius/P.radius)^2 +
                                                a6 * (S.radius/P.radius)^3 +
                                                a7 * (S.radius/P.radius)^4))*
          exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-sigma*pi*P.radius^2)/kbT)*pi/kbT}      
      
      Hp5<-function(P.radius) {P.radius^5*(6 * pi * (1-(S.radius/P.radius))^2 /
                                             (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                                (1 + a1 * (1-(S.radius/P.radius)) + 
                                                   a2 * (1-(S.radius/P.radius))^2) +
                                                a3 +
                                                a4 * (S.radius/P.radius) +
                                                a5 * (S.radius/P.radius)^2 +
                                                a6 * (S.radius/P.radius)^3 +
                                                a7 * (S.radius/P.radius)^4))*
          exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-sigma*pi*P.radius^2)/kbT)*(pi/kbT)*(2*pi/kbT)}
      
      Hp6<-function(P.radius) {P.radius^6*(6 * pi * (1-(S.radius/P.radius))^2 /
                                             (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                                (1 + a1 * (1-(S.radius/P.radius)) + 
                                                   a2 * (1-(S.radius/P.radius))^2) +
                                                a3 +
                                                a4 * (S.radius/P.radius) +
                                                a5 * (S.radius/P.radius)^2 +
                                                a6 * (S.radius/P.radius)^3 +
                                                a7 * (S.radius/P.radius)^4))*
          exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-sigma*pi*P.radius^2)/kbT)*(pi/kbT)^2}
      I2<-integrate(Hp2, lower = S.radius, upper = gamma/sigma)
      I3<-integrate(Hp3, lower = S.radius, upper = gamma/sigma)
      I4<-integrate(Hp4, lower = S.radius, upper = gamma/sigma)
      I5<-integrate(Hp5, lower = S.radius, upper = gamma/sigma)
      I6<-integrate(Hp6, lower = S.radius, upper = gamma/sigma)
      Imatrix[1,j]<-I2$value
      Imatrix[2,j]<-I3$value
      Imatrix[3,j]<-I4$value
      Imatrix[4,j]<-I5$value
      Imatrix[5,j]<-I6$value
    }
    #
    #
    # perhaps this parameter is the key to solve the problem
    if (iter ==0 & plot.flag){scale <-min(abs(DTD.scaled[1])/Imatrix[1,]) }#  else{
#      scale <- scale+alpha*(min(abs(DTD.scaled)/Imatrix[1,])-scale)
#    }# this scaling factor may be the problem...
    
    
    I2<-scale*Imatrix[1,]
    I3<-scale*Imatrix[2,]
    I4<-scale*Imatrix[3,]
    I5<-scale*Imatrix[4,]
    I6<-scale*Imatrix[5,]

    A <- matrix(c(sum(I3*I3*kbT+4*pi*(I2-DTD.scaled)*I4),-sum(I3*I4*kbT+(I2-DTD.scaled)*I5*kbT),
                  sum(I3*I4*kbT+(I2-DTD.scaled)*I5*kbT),-sum(I4*I4*kbT+(I2-DTD.scaled)*I6*kbT)),nrow = 2,byrow=TRUE)
#    A <- -sum(I4*I4*kbT+(I2-DTD.scaled)*I6*kbT)
    invA <- matrix(c(-sum(I4*I4*kbT+(I2-DTD.scaled)*I6*kbT),sum(I3*I4*kbT+(I2-DTD.scaled)*I5*kbT),
                     -sum(I3*I4*kbT+(I2-DTD.scaled)*I5*kbT),sum(I3*I3*kbT+4*pi*(I2-DTD.scaled)*I4)),nrow = 2,byrow=TRUE)/det(A)
    
    product <- c(sum((I2-DTD.scaled)*I3*kbT),sum((I2-DTD.scaled)*I4*kbT))
#    product <- sum((I2-DTD.scaled)*I4*kbT)
    
    delta <-invA %*% product
#    delta <- product/A
#    
#    gamma<-gamma+alpha*delta[1]
#    sigma<-min(sigma+alpha*delta,gamma/max(S.radii)*0.999)
    sigma<- sigma +alpha*delta[2]
    #gamma<- max(sigma*max(S.radii)/0.999, gamma+alpha*delta[1])
    #print(c(delta[1]/gamma,",",delta[2]/sigma) )
    #eps <- max(abs(delta[1]/gamma),abs(delta[2]/sigma))
    #eps <- max(abs(delta[1]/gamma),abs(delta[2]/sigma))
    eps <- abs(delta[2]/sigma)
    if(mod(iter,iter_max/100) ==0){
      DTD.data <- data.frame(t(rbind(S.radii,DTD.scaled,I2)))
      colnames(DTD.data)<-c("S.radii","DTD.scaled","DTD.estimated")
      DTD.data.ggplot <-melt(DTD.data,id.vars  = "S.radii")
      print(ggplot(DTD.data.ggplot,aes(x=S.radii,y=value,color=variable))+geom_point())
      print(c(iter, abs(delta[1]/gamma),abs(delta[2]/sigma), gamma/sigma,sigma, gamma, sum((I2-DTD.scaled))))
      #print(c(iter, abs(delta/sigma), gamma/sigma,sigma, gamma))
      delta_history[100*iter/iter_max+1,]<-c(iter,delta[1]/gamma,delta[2]/sigma)
      par_history[100*iter/iter_max+1,]<-c(iter,gamma,sigma,gamma/sigma)
    }
    iter <- iter+1
  }
  DTD.data <- data.frame(t(rbind(S.radii,DTD.scaled,I2)))
  colnames(DTD.data)<-c("S.radii","DTD.scaled","DTD.estimated")
  DTD.data.ggplot <-melt(DTD.data,id.vars  = "S.radii")
  ggplot(DTD.data.ggplot,aes(x=S.radii,y=value,color=variable))+geom_point()
  delta_history<-data.frame(delta_history)
  colnames(delta_history)<-c("iter","dgamma","dsigma")
  par_history<-data.frame(par_history)
  colnames(par_history)<-c("iter","gamma","sigma","radius")
  
  if (plot.flag){
  print(
    ggplot(par_history,aes(x=iter,y=gamma))+geom_point()+
  ggplot(melt(delta_history,id.vars="iter"),aes(x=iter,y=abs(value),color=variable))+geom_point()+scale_y_log10()+
    ggplot(par_history,aes(x=iter,y=sigma))+geom_point()+
  ggplot(par_history,aes(x=iter,y=radius))+geom_point())
  }
  return(c(gamma,sigma,scale))
}
#
# compute multiple radii
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
        weighted.count <- DTD.scaled[,t]/fm.predict
        print(weighted.count)
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

