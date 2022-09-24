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
# FreemanSA,WangMA,WeaverJC.1994.Theory of electroporation for a planar bilayer membrane:
# predictions of the fractional aqueous area, change in capacitance and pore-pore separation.
# Biophys. J. 67:42–56
gamma<-2e-12
# radius of BSA
S.radius<- 3.48e-9
#
sigma.list <- seq(10,200,length.out = 50)*1e-6
# data matrix
I2<-matrix(nrow=length(sigma.list),ncol=1)

# compute flux
for (icnt in 1:length(sigma.list)){
  sigma<-sigma.list[icnt]
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
  
  I2_1 <-integrate(Hp2, lower = S.radius, upper = gamma/sigma)
  I2[icnt,1]<-I2_1$value
}

I2<-data.frame(I2)
I2$sigma<-sigma.list
#I2$radius <- gamma/I2$sigma
I2$type<-"simulation"
colnames(I2)<-c("amount","surface_tension","type")

ggplot(I2,aes(x=surface_tension,y=amount))+geom_line()+theme_bw()

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


