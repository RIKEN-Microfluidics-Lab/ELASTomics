elp <-read.table(file = "/home/samba/public/shintaku/github/ELASTomics/DTD_elp.csv",sep = ",",header = TRUE)
stokes <-read.table(file = "/home/samba/public/shintaku/github/ELASTomics/DTD_stokes.csv",sep = ",",header = TRUE)

ggplot(elp,aes(x=MW,y=mean,shape=type,color=type))+geom_point(size=2)+scale_x_log10()+
  geom_errorbar(aes(ymin = mean - 0.5*standard_error, ymax = mean + 0.5*standard_error),width=0.05) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Molecular weight of dextran (kDa)",y=bquote('Electrophoretic mobility'~(m^2/Vs)),size=20)+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12))+
ggplot(stokes,aes(x=MW,y=0.5*mean,shape=type,color=type))+geom_point(size=2)+scale_x_log10(limits= c(3,600))+
  geom_errorbar(aes(ymin = 0.5*mean - 0.25*standard_error, ymax = 0.5*mean + 0.25*standard_error),width=0.05) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Molecular weight of dextran (kDa)",y="Stokes radius (nm)",size=20)+
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12))+
  scale_x_log10(limits= c(3,600))+ylim(c(0,15))
        