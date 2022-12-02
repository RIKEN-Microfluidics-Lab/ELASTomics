DimPlot(cancer_merge1234,reduction="tsne")
Idents(cancer_merge1234)<-cancer_merge1234[["celltype"]]
FeaturePlot(subset(cancer_merge1234,subset=nCount_RNA>3000 & percent.mt<10),features="percent.mt", reduction = "tsne")

#
# subset useful cell data and visualize
cancer_subset <- subset(cancer_merge1234,subset=nCount_RNA>3000 & percent.mt<10 & condition=="normal")
cancer_subset<-subset(cancer_subset,subset=NEP=="75V",invert=TRUE)
cancer_subset <- FindVariableFeatures(cancer_subset, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_subset)
cancer_subset <- ScaleData(cancer_subset, features = all.genes)
cancer_subset <- RunPCA(cancer_subset, npcs=20, features = VariableFeatures(object = cancer_subset))
#DimPlot(cancer_subset,group.by = "celltype")
cancer_subset <- RunUMAP(cancer_subset, dims = 1:10)
cancer_subset <- RunTSNE(cancer_subset,dims=1:10)

DimPlot(cancer_subset,reduction="tsne")+
  DimPlot(cancer_subset,reduction = "tsne",group.by = "NEP")+
  FeaturePlot(cancer_subset,reduction = "tsne",features = "percent.mt")

p1<-RidgePlot(subset(cancer_subset,subset=celltype=="MCF10A"),
              features = "dtd_FLD004",group.by = "NEP",cols=c("#D9D7D7","#F8766D"))+
  labs(title = "MCF10A")+xlim(c(0,8))
p2<-RidgePlot(subset(cancer_subset,subset=celltype=="MCF7"),
              features = "dtd_FLD004",group.by = "NEP",cols=c("#D9D7D7","#7CAE00"))+
  labs(title = "MCF7")+xlim(c(0,8))
p3<-RidgePlot(subset(cancer_subset,subset=celltype=="MDAMB231"),
              features = "dtd_FLD004",group.by = "NEP",cols=c("#D9D7D7","#00BFC4"))+
  labs(title = "MDA-MB-231")+xlim(c(0,8))
p4<-RidgePlot(subset(cancer_subset,subset=celltype=="PC3"),
              features = "dtd_FLD004",group.by = "NEP",cols=c("#D9D7D7","#C77CFF"))+
  labs(title = "PC-3")+xlim(c(0,8))
gridExtra::grid.arrange(p1,p2,p3,p4,nrow=1)

mcf10a <- subset(cancer_subset, celltype=="MCF10A")
FeaturePlot(mcf10a,features="dtd_FLD004",reduction="tsne")+DimPlot(mcf10a,reduction="tsne",group.by="NEP")
RidgePlot(mcf10a,features = c("dtd_FLD004","dtd_FLD010","dtd_FLD040","dtd_FLD070","dtd_FLD150","dtd_FLD500"),
          group.by = "NEP")

fldmean<-c(AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD004",slot="counts")$DTD,
  AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD010",slot="counts")$DTD,
  AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD040",slot="counts")$DTD,
  AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD070",slot="counts")$DTD,
  AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD150",slot="counts")$DTD,
  AverageExpression(subset(mcf10a,subset=NEP=="0V"),features="FLD500",slot="counts")$DTD)

fldcount <- data.frame(mcf10a[["DTD"]]@counts)
mcf10a_dtd<-data.frame(fldcount[1:6,])
fldcount <- t(sweep(fldcount[1:6,],MARGIN=2,fldmean,FUN = "-"))
library(reshape2)
ggplot(melt(fldcount),aes(x=Var2,y=value-min(fldcount)))+geom_violin()+geom_jitter(size=0.1)+scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y="Number of transported DTD (counts)",x="Molecular weight of dextran (kDa)",size=20)
fldcount<-t(fldcount)
mcf10a_bk <- CreateAssayObject(sweep(fldcount[1:6,],MARGIN=2,fldmean,FUN = "-")-min(fldcount))
mcf10a[["DTDbk"]]<-mcf10a_bk
RidgePlot(mcf10a,features = "dtdbk_FLD004",group.by = "NEP")+scale_x_log10()
RidgePlot(mcf10a,features = c("dtdbk_FLD004","dtdbk_FLD010",
                              "dtdbk_FLD040","dtdbk_FLD070",
                              "dtdbk_FLD150","dtdbk_FLD500"),
          group.by = "NEP")+scale_x_log10()
par<-data.frame(t(c(2e-11,5.359301e-06)))
colnames(par)<-c("gamma","sigma")


output <- elast.comp.radius(mcf10a_dtd$cancer4TTCGCTGCAATGCAGG,S.radii,par,scale,plot.flag=TRUE)


