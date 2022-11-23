library(pracma)
library(stats)
library(sys)
library(matlib)
library(numbers)
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
#KotnikT,RemsL,TarekM,MiklavcicD,
#Membrane elecgtroporation and electropermeabilization: mechanism and models
#Annual Review of Biophysics
d<-5e-9
epse<-7.1e-10
epsm<-4.4e-11

# FreemanSA,WangMA,WeaverJC.1994.Theory of electroporation for a planar bilayer membrane:
# predictions of the fractional aqueous area, change in capacitance and pore-pore separation.
# Biophys. J. 67:42–56
gamma<-2e-11
# radius of BSA
S.radius<- 3.48e-9
#
sigma<-1e-4
V<-50.0e-3
sigma.list <- seq(10,800,length.out = 5)*1e-6
V.list <- seq(10,200,length.out = 5)*1e-3
P.radius <-seq(0.4,20,length.out=100)*1e-9
#
# Boltzmann distribution with various voltages/surface tension
for (icnt in 1:length(sigma.list)){
  #V<-V.list[icnt]
  sigma<-sigma.list[icnt]
  z<-function(radius){exp(-(2*pi*gamma*radius*(1+C/radius^5)-(sigma+V^2*(epse-epsm)/2/d)*pi*radius^2)/kbT)}
  Z<-integrate(z, lower = 0, upper = max(P.radius))
  prob<-exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-(sigma+V^2*(epse-epsm)/2/d)*pi*P.radius^2)/kbT)/(Z$value/max(P.radius))
  prob<-data.frame(prob)
  prob$radius <-P.radius
  prob$sigma<-sigma
  prob$voltage <-V
  prob$free_energy<-(2*pi*gamma*P.radius*(1+C/P.radius^5)-(sigma+V^2*(epse-epsm)/2/d)*pi*P.radius^2)/kbT
  if (icnt==1){prob_stack<-prob}
  else
  {
      prob_stack<-rbind(prob_stack,prob)
    }
}
p2<-ggplot(prob_stack,aes(x=radius*1e9,y=prob,color=as.factor(sigma)))+geom_line()+scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1<-ggplot(prob_stack,aes(x=radius*1e9,y=free_energy,color=as.factor(sigma)))+geom_line()+#scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#
# amount of imported molecules with the radius of S.radius
#
S.list<-seq(3,18,length.out = 20)*1e-9
V.list<-seq(50,50,length.out = 20)*1e-3
I2<-matrix(nrow=length(V.list),ncol=1)
# compute flux
for (jcnt in 1:length(sigma.list)){
  sigma<-sigma.list[jcnt]
  for (icnt in 1:length(S.list)){
    V<-V.list[icnt]
    S.radius<-S.list[icnt]
    Hp2<-function(P.radius) {P.radius^2*(6 * pi * (1-(S.radius/P.radius))^2 /
                                           (b1 * (1-(S.radius/P.radius))^(-5/2) *
                                              (1 + a1 * (1-(S.radius/P.radius)) +
                                                 a2 * (1-(S.radius/P.radius))^2) +
                                              a3 +
                                              a4 * (S.radius/P.radius) +
                                              a5 * (S.radius/P.radius)^2 +
                                              a6 * (S.radius/P.radius)^3 +
                                              a7 * (S.radius/P.radius)^4))*
        exp(-(2*pi*gamma*P.radius*(1+C/P.radius^5)-(sigma+V^2*(epse-epsm)/2/d)*pi*P.radius^2)/kbT)}
    
    I2_1 <-integrate(Hp2, lower = S.radius, upper = gamma/(sigma+V^2*(epse-epsm)/2/d))
    I2[icnt,1]<-I2_1$value*V
  }
  
  I2<-data.frame(I2)
  I2$sigma<-sigma
  I2$V <- V.list
  I2$r <- S.list
  #I2$radius <- gamma/I2$sigma
  I2$type<-"simulation"
  if(jcnt==1){I2_stack<-I2}
  else{
    I2_stack<-rbind(I2_stack,I2)
  }
}

#colnames(I2)<-c("amount","surface_tension","type")
colnames(I2_stack)<-c("amount","surface_tension","voltage","radius","type")

ggplot(I2_stack,aes(x=radius*1e9,y=amount,color=as.factor(surface_tension)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line()+scale_y_log10()#+
#  coord_cartesian(xlim = NULL, ylim =c(1e-100,1e-68))
p1+p2+p3





exp.data <-read.table("/home/samba/public/shintaku/ELASTomics/SD2_500nm_ResultData.csv",sep = ",",header = TRUE)
base<-min(exp.data$RawIntDen)
scale<-max(exp.data$RawIntDen-base)/max(I2$amount)
exp.data.sub<-data.frame(cbind((exp.data$RawIntDen-base)/scale,exp.data$Med_ST*1e-6))
exp.data.sub$type <-"exp"
colnames(exp.data.sub)<-c("amount","surface_tension","type")

data<-rbind(I2,exp.data.sub)
cor_val <-as.character(cor(exp.data.sub$amount,exp.data.sub$surface_tension,method = "spearman"))
ggplot()+geom_point(data=exp.data.sub,aes(x=surface_tension,y=amount,color=type))+theme_bw()+
  geom_line(data=I2,aes(x=surface_tension,y=amount,color=type))+scale_x_log10()+ylim(c(0,1e-32))

