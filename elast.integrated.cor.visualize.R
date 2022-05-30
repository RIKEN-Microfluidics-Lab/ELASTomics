
library(Seurat)
library(ggplot2)
p1<-RidgePlot(tig.combined.nep,features = "dtd_FLD004",group.by = "age")+scale_fill_manual(values = c( "#0072B2","#D55E00"))
p2<-RidgePlot(tig.combined.nep,features = "dtd_FLD500",group.by = "age")+ scale_fill_manual(values = c( "#0072B2","#D55E00"))
p3<-RidgePlot(tig.combined,features = "dtd_FLD004",group.by = "condition")+scale_fill_manual(values = c( "#0072B2","#D55E00"))
p1+p2+p3

#aging related selected from knowledge
#
p1<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "ATF3",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p2<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "FOS",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p3<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "ADM",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p4<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "CDKN1A",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p5<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "DTYMK",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p6<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "NCL",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p7<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "PCLAF",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p8<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "DLGAP5",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
anti_aging_genes <-c("DTYMK","NCL","PCLAF","DLGAP5")
gridExtra::grid.arrange(p1, p2,p3, p4, p5, p6,p7, p8, nrow = 2)

#
# RRAD and RHOB
p1<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "RRAD",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p2<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "RHOB",pt.size = 0.1,group.by = "age")+scale_color_manual(values = c( "#0072B2","#D55E00"))
p1+p2

#
#https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13315 
#
FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IGFBP6",pt.size = 0.1,group.by = "age")+
  FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "KLF2",pt.size = 0.1,group.by = "age")+
  FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "MYBL2",pt.size = 0.1,group.by = "age")
#SASP atlas genes
p1<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "TIMP1",pt.size = 0.1,group.by = "age")
p2<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "SERPINE2",pt.size = 0.1,group.by = "age")
p3<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "RARRES2",pt.size = 0.1,group.by = "age")
p4<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "PLAT",pt.size = 0.1,group.by = "age")
p5<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IGFBP6",pt.size = 0.1,group.by = "age")
p6<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IGFBP5",pt.size = 0.1,group.by = "age")
p7<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IGFBP4",pt.size = 0.1,group.by = "age")
p8<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IGFBP2",pt.size = 0.1,group.by = "age")
p9<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "GDF15",pt.size = 0.1,group.by = "age")
p10<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "FN1",pt.size = 0.1,group.by = "age")
p11<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "F3",pt.size = 0.1,group.by = "age")
p12<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "DCN",pt.size = 0.1,group.by = "age")
p13<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "CTHRC1",pt.size = 0.1,group.by = "age")
p14<-FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "COL3A1",pt.size = 0.1,group.by = "age")
gridExtra::grid.arrange(p1, p2,p3, p4, p5, p6,p7, p8, p9, p10,p11,p12,p13,p14,nrow = 3)







