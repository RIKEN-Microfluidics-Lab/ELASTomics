library(ggplot2)
library(ggrepel)
library(Seurat)
library(pheatmap)
library(destiny)
#DoHeatmap(tig.bulk.seurat,features=aging_gene)
exp.matrix<-data.frame(tig.combined.ctl[["RNA"]]@data)
age <- data.frame(tig.combined.ctl[["age"]])
rownames(age)<-gsub("-",".",rownames(age))

exp.matrix<-exp.matrix[rownames(exp.matrix) %in% aging_gene,]
pheatmap(exp.matrix,scale="row")
dm <- DiffusionMap(data.frame(t(exp.matrix)))

pheatmap(exp.matrix[,order(dm$DC1,decreasing = FALSE)],
         scale="row",
         cluster_cols = FALSE,
         annotation_col = age,
         show_colnames = FALSE)
age$DC1<-dm$DC1
ggplot(age,aes(x=age,y=-DC1))+geom_violin()+geom_jitter(size=0.1)

#tig.combined.nep <- RenameIdents(tig.combined.nep, `0` = "TIG1-50", `1` = "TIG1-20")
#
# 
FeaturePlot(tig.combined.nep,features=c("dtd_FLD004","ATF3","TNFAIP3","CXCL8"),reduction="pca")
FeaturePlot(tig.combined.ctl,features=c("dtd_FLD004","ATF3","TNFAIP3","CXCL8"),reduction="pca")
VlnPlot(tig.combined,features = c("ATF3","TNFAIP3","CXCL1","CXCL2","CXCL3","CXCL8","FOS"),
        group.by = "condition",slot="data",pt.size=0)

# Prof. Akiko Takahashi senescence TIG-3
# https://www.nature.com/articles/s41467-018-03555-8#Fig1
VlnPlot(tig.combined.ctl,features = c("FOS","CDKN2A","CGAS"),group.by = "age",slot="data",pt.size=1)
VlnPlot(tig.combined.nep,features = c("FOS","CDKN2A","CGAS"),group.by = "age",slot="data",pt.size=1)
# virus infection DNase related
VlnPlot(tig.combined,features = c("CGAS","IFNB1","AIM2","AMPD1","STING1","FOS"),group.by = "condition",slot="data",pt.size=1)

DimPlot(tig.combined.nep)+FeaturePlot(tig.combined.nep,features=c("dtd_FLD040"))
FeaturePlot(tig.combined.nep,features=c("SAT1","FOS","FOSB","ADM"))
FeaturePlot(tig.combined.nep,features=c("TMPO","DTYMK","NCL","PCLAF"))


#source("elast.CellCycle.TIG1.R")