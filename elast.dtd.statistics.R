library(reshape2)
S.radii <- c(1.4, 2.7, 6.3, 15.1)
tig.nep <- NormalizeData(tig.nep,assay="DTD",normalization.method = "LogNormalize",scale.factor = 1)
VlnPlot(tig.nep,features=c("dtd_FLD004","dtd_FLD010","dtd_FLD040","dtd_FLD070","dtd_FLD150","dtd_FLD500"),slot="counts")
VlnPlot(tig.nep,features=c("dtd_FLD004","dtd_FLD010","dtd_FLD040","dtd_FLD070","dtd_FLD150","dtd_FLD500"))
DimPlot(tig.nep)+FeaturePlot(tig.nep,features="dtd_FLD004",slot="data")+FeaturePlot(tig.nep,features="dtd_FLD500",slot="data")

tig.dtd.scale <-tig.nep[["DTD"]]@counts
tig.dtd.scale<-sweep(as.matrix(tig.dtd.scale),1,as.matrix(concentration),"/")
tig.dtd.scale <-  t(tig.dtd.scale[1:6,])
tig.dtd.scale.melt <- melt(tig.dtd.scale)
colnames(tig.dtd.scale.melt) <-  c("cellid","DTD","Scaled.count")
ggplot(tig.dtd.scale.melt,aes(x=DTD,y=Scaled.count))+geom_violin()+geom_jitter(size=0.1)+scale_y_log10()

colnames(Res)<- c("radii","se","t","pval")
ggplot(Res,aes(x=radii))+ geom_histogram()+scale_x_log10()
ggplot(Res,aes(x=pval))+ geom_histogram()+scale_x_log10()

tig.50<-subset(tig.nep,idents = "Unknown")
median(as.numeric(unlist(tig.50[["radii"]])))
