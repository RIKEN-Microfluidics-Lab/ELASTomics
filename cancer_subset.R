cancer_normal <-subset(subset(cancer_merge123456,subset=condition=="normal"),subset=NEP=="75V",invert=TRUE)
DimPlot(cancer_normal,reduction="tsne",group.by = "celltype")+DimPlot(cancer_normal,reduction="tsne",group.by = "NEP")+
  DimPlot(cancer_normal,reduction="tsne",group.by = "run")

Idents(cancer_normal)<-cancer_normal[["celltype"]]
FeaturePlot(cancer_normal,features="nCount_RNA", reduction = "tsne")

#
# subset useful cell data and visualize
cancer_subset <- subset(cancer_normal,subset=nCount_RNA>3000 & percent.mt<5 & condition=="normal")
#cancer_subset<-subset(cancer_subset,subset=NEP=="75V",invert=TRUE)
cancer_subset <- FindVariableFeatures(cancer_subset, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_subset)
cancer_subset <- ScaleData(cancer_subset, features = all.genes)
cancer_subset <- RunPCA(cancer_subset, npcs=20, features = VariableFeatures(object = cancer_subset))
#DimPlot(cancer_subset,group.by = "celltype")
cancer_subset <- RunUMAP(cancer_subset, dims = 1:10)
cancer_subset <- RunTSNE(cancer_subset,dims=1:10)

DimPlot(cancer_subset,reduction="umap")+
  DimPlot(cancer_subset,reduction = "umap",group.by = "NEP")+
  FeaturePlot(cancer_subset,reduction = "umap",features = "percent.mt")

p1<-RidgePlot(subset(cancer_subset,subset=celltype=="MCF10A"),
              features = "dtd_FLD004",group.by = "NEP",cols=c("#D9D7D7","#F8766D"))+
  labs(title = "MCF10A")+xlim(c(0,8))
RidgePlot(subset(cancer_subset,subset=celltype=="MCF7" & run==c("fifth","sixth")),
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

mcf10aNEP <- subset(mcf10a,subset = NEP=="40V")

DTDcoc <- c(10455, 9811, 10015, 4647, 11254, 6004)
DTDmob <- c(2.26E-08,2.29E-08,2.38E-08,1.72E-08,1.44E-08,1.69E-08)
DTDcoc*DTDmob

fldcount <- t(data.frame(mcf10aNEP[["DTD"]]@counts))
#mcf10a_dtd<-data.frame(fldcount[1:6,])
#fldcount <- t(sweep(fldcount[1:6,],MARGIN=2,fldmean,FUN = "-"))
library(reshape2)
ggplot(melt(fldcount[,1:6]*DTDcoc*DTDmob),aes(x=Var2,y=value-min(fldcount)))+geom_violin()+geom_jitter(size=0.1)+scale_y_log10()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y="Number of transported DTD (counts)",x="Molecular weight of dextran (kDa)",size=20)

fldmatrix <- data.frame(fldcount[,1:6]*DTDcoc*DTDmob)
fldmatrix$ratio <- fldmatrix$FLD040/fldmatrix$FLD004
ggplot(fldmatrix,aes(x=FLD004,y=FLD040))+geom_point()
cor(fldmatrix$FLD004,fldmatrix$FLD500)

fldcount<-t(fldcount)
mcf10a_bk <- CreateAssayObject(sweep(fldcount[1:6,],MARGIN=2,fldmean,FUN = "-")-min(fldcount))
mcf10a[["DTDbk"]]<-mcf10a_bk
RidgePlot(mcf10a,features = "dtdbk_FLD004",group.by = "NEP")+scale_x_log10()
RidgePlot(mcf10a,features = c("dtdbk_FLD004","dtdbk_FLD010",
                              "dtdbk_FLD040","dtdbk_FLD070",
                              "dtdbk_FLD150","dtdbk_FLD500"),
          group.by = "NEP")+scale_x_log10()
#par<-data.frame(t(c(2e-11,1e-02)))
par<-data.frame(t(c(1.138212e-13,4.859301e-06)))
colnames(par)<-c("gamma","sigma")


output <- elast.comp.radius(mcf10a_dtd[seq(1,4),]$cancer4AAGCATCAGCCTAACT,S.radii[seq(1,4)],par,1e-6,scale,plot.flag=TRUE)

mcf10a_dtd$cancer4AAGCATCAGCCTAACT
