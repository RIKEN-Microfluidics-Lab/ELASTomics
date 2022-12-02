require(readr)
library(plyr,dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse,R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(Seurat)
library(SingleCellSignalR)
library(seqinr)
library(stringr)
refexp <- read.csv("/home/samba/sanger/shintaku/ELASTomics/Expression_22Q2_Public_subsetted.csv")
refexp<-data.frame(t(refexp))
celltype<-refexp[2,]
expression<-data.frame(as.matrix(refexp[7:nrow(refexp),]))
colnames(expression)<-celltype
expression$MCF7_PC3 <- as.numeric(expression$MCF7)-as.numeric(expression$PC3)

rdir <- "/home/samba/public/shintaku/github/ELASTomics/"
#
# load ELASTomics data from an output of cellranger with cite-seq pipeline
#
source("elast.load.elast.data.R")
# sample data sheet
# HiSeqX_010_011
# https://docs.google.com/spreadsheets/d/1GJjzlQfBgJoL_iBLtggQhTF5ABsY_ZZz/edit?rtpof=true#gid=1765377542
# HiSeqX_012_013
# https://docs.google.com/spreadsheets/d/1vO2SCCNvDIM0N9YOOMDDsC6Xiqjavxqi/edit#gid=1765377542

#First
# ・PC3 control
# ・PC3 EP(40V) ELPOligoNo.35
# ・MCF7 EP (40V) ELPOligoNo.37
#Second
# ・PC3 CytoD EP(40V) ELPOligoNo.36
# ・MCF7 Control
#Third
# ・MDAMB231 EP(40V) ELPOligoNo.32
# ・PC3 EP(75V) ELPOligoNo.33
# ・MCF10A Control
#Forth
# ・MCF10A EP(40V) ELPOligoNo.31
# ・PC3 MBCD EP(40V) ELPOligoNo.34
# ・MDAMB231 Control
#Fifth
# ・MCF10A EP(40V) ELPOligoNo.31
# ・MCF7 Control
#Sixth
# ・MCF7 EP (40V) ELPOligoNo.37
# ・MCF10A CytoD EP(40V) ELPOligoNo.36

wdir <- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX010_011/cancer1_12/outs/filtered_feature_bc_matrix/"
cancer1 <-load.elast.data(wdir,"cancer1",100)
cancer1[["run"]]<-"first"
wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX010_011/cancer2_12/outs/filtered_feature_bc_matrix/"
cancer2 <-load.elast.data(wdir,"cancer2",100)
cancer2[["run"]]<-"second"
cancer <-cancer2
source("cancer2_Seurat.clustering.R")
cancer <- merge(cancer1,y=cancer2)
source("cancer1_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX012_013/cancer3_12/outs/filtered_feature_bc_matrix/"
cancer3 <-load.elast.data(wdir,"cancer3",100)
cancer3[["run"]]<-"third"
cancer<-merge(cancer3,y=cancer_merge)
source("cancer3_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX012_013/cancer4_12/outs/filtered_feature_bc_matrix/"
cancer4 <-load.elast.data(wdir,"cancer4",100)
cancer4[["run"]]<-"forth"
cancer<-merge(cancer4,y=cancer_merge123)
source("cancer4_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20221012HiSeqX016/MCF7_10A/outs/filtered_feature_bc_matrix/"
cancer5 <-load.elast.data(wdir,"cancer5",100)
cancer5[["run"]]<-"fifth"
cancer<-cancer5
source("cancer5_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20221012HiSeqX017/MCF7_10A/outs/filtered_feature_bc_matrix/"
cancer6 <-load.elast.data(wdir,"cancer6",100)
cancer6[["run"]]<-"sixth"
cancer<-cancer6
source("cancer6_Seurat.clustering.R")

VlnPlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF10A"),features = "dtd_FLD004",group.by = "condition", cols = c("#D55E00", "#0072B2"))
VlnPlot(subset(cancer_merge123456,subset=celltype=="MCF7"),features = "dtd_FLD004",group.by = "NEP", cols = c("#D55E00", "#0072B2"))
VlnPlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth"), invert=TRUE),subset=celltype=="PC3"),features = "dtd_FLD004",group.by = "NEP", cols = c("#D55E00","#40A39A", "#0072B2")) 
RidgePlot(subset(subset(cancer_merge1234,subset=NEP=="40V"),subset=condition=="normal"), features = "dtd_FLD004",group.by = "celltype")
RidgePlot(subset(subset(cancer_merge123456,subset=run==c("first", "second")),subset=celltype=="MCF7"),features = "dtd_FLD004",group.by = "NEP")+
  RidgePlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF7"),features = "dtd_FLD004",group.by = "NEP")
DimPlot(subset(subset(cancer_merge123456,subset=run==c("first", "second")),subset=celltype=="MCF7") ,reduction="pca" ,group.by = "NEP")+
  DimPlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF7") ,reduction="pca" ,group.by = "NEP")

# Dextran ratio
DTDcoc <- c(10455, 9811, 10015, 4647, 11254, 6004)
DTDmob <- c(2.26E-08,2.29E-08,2.38E-08,1.72E-08,1.44E-08,1.69E-08)
DTD40V <- subset(subset(cancer_merge123456,subset=NEP=="40V"),subset=condition=="normal")
DTD40V_dtd <-data.frame(t(DTD40V[["DTD"]]@counts))
DTD40V_dtd$FLD004 <- DTD40V_dtd$FLD004 * (10000 / DTDcoc[1]) / (DTDmob[1] / 1E-07)
DTD40V_dtd$FLD010 <- DTD40V_dtd$FLD010 * (10000 / DTDcoc[2]) / (DTDmob[2] / 1E-07)
DTD40V_dtd$FLD040 <- DTD40V_dtd$FLD040 * (10000 / DTDcoc[3]) / (DTDmob[3] / 1E-07)
DTD40V_dtd$FLD070 <- DTD40V_dtd$FLD070 * (10000 / DTDcoc[4]) / (DTDmob[4] / 1E-07)
DTD40V_dtd$FLD150 <- DTD40V_dtd$FLD150 * (10000 / DTDcoc[5]) / (DTDmob[5] / 1E-07)
DTD40V_dtd$FLD500 <- DTD40V_dtd$FLD500 * (10000 / DTDcoc[6]) / (DTDmob[6] / 1E-07)
DTD40V_dtd_melt<-melt(DTD40V_dtd,measure.vars = c("FLD004","FLD010","FLD070","FLD500"),value.name = "Counts")
ggplot(DTD40V_dtd_melt,aes(x=variable,y=Counts))+geom_violin()+geom_jitter(size=0.1)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")

source("cancer.integraated.R")

DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "celltype")
DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "NEP", cols = c("#D55E00", "#0072B2"))
DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "run")
FeaturePlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,features = "FLD004", max.cutoff = 4)
FeaturePlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,features = "FLD500", max.cutoff = 4)
VlnPlot(subset(cancer.integrated,subset=NEP=="40V"), features = "FLD004", group.by = "celltype")
VlnPlot(subset(cancer.integrated,subset=NEP=="40V"), features = "FLD500", group.by = "celltype")

#DNAtag("dtd_M1AE31", "dtd_MM2E32", "dtd_P3HE33", "dtd_P3ME34", "dtd_P3NE35", "dtd_P3CE36", "dtd_MC7E37")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_M1AE31",reduction = "tsne", max.cutoff = 6) #MCF10Atag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_MM2E32",reduction = "tsne", max.cutoff = 2) #MDAMB231tag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_P3NE35",reduction = "tsne", max.cutoff = 2, min.cutoff = 0.5) #PC3tag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_MC7E37",reduction = "tsne", max.cutoff = 4, min.cutoff = 1) #MCF7tag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"),
            features = c("dtd_M1AE31", "dtd_MM2E32", "dtd_P3HE33", "dtd_P3ME34", "dtd_P3NE35", "dtd_P3CE36", "dtd_MC7E37"),
            reduction = "tsne", max.cutoff = 4, min.cutoff = 2)

FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "tsne", features = "CYBA", max.cutoff = 6) #MCF7
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "tsne", features = "S100A2", max.cutoff = 6) #PC-3
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "tsne", features = "CAVIN3", max.cutoff = 4) #MDAMB231
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "tsne", features = "DUSP1", max.cutoff = 4) #MCF10A

#FindMarkers(subset(cancer.integrated,subset=NEP=="40V"), ident.1 = "MCF10A")
#
# cell cycle
#
gene_list <- data.frame(unique(rownames(cancer.integrated)))
source(file.path(rdir,"elast.biomaRt.R"))
source(file.path(rdir,"util/fucci_cellcycle_genes.R"))
  
sub_ref <- ref %>%
  dplyr::filter(gene_short_name %in% rownames(cancer.integrated))
  
genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
VlnPlot(cancer.integrated,features = "dtd_FLD004")+
  FeatureScatter(cancer.integrated,feature2 = "dtd_FLD004",feature1 = "S.Score")+
  FeatureScatter(cancer.integrated,feature2 = "dtd_FLD004",feature1 = "G2M.Score")+
  FeatureScatter(cancer.integrated,feature1 = "G2M.Score" ,feature2 = "S.Score")

# mcf10a subset 
mcf10a <- subset(cancer.integrated,subset=celltype=="MCF10A")
Idents(mcf10a)<-mcf10a[["NEP"]]
mcf10a<-NormalizeData(mcf10a,assay="RNA",normalization.method = "LogNormalize", scale.factor = 1e5)
mcf_nep <- FindMarkers(mcf10a,ident.1 = "40V",ident.2 = "0V")
View(mcf_nep)

gene.list<-symbol2entrez(rownames(mcf_nep[mcf_nep$avg_log2FC>1 & mcf_nep$p_val_adj<1e-3,]))

ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.01, 
                       readable      = TRUE)
barplot(ego_result)
View(data.frame(ego_result))

gene.plus.list<-symbol2entrez(rownames(mcf_nep[mcf_nep$avg_log2FC>0.5 & mcf_nep$p_val_adj<1e-3,]))
gene.minus.list<-symbol2entrez(rownames(mcf_nep[mcf_nep$avg_log2FC< -0.5 & mcf_nep$p_val_adj<1e-3,]))
gene.list <- list(negative=gene.minus.list$gene_id,
                  positive=gene.plus.list$gene_id)
xx <- compareCluster(gene.list, fun="enrichKEGG")
summary(xx)
dotplot(xx)
#

mcf10a.list <- SplitObject(cancer_merge, split.by = "NEP")
mcf10a.list <- lapply(X =  mcf10a.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = mcf10a.list)

mcf10a.anchors <- FindIntegrationAnchors(object.list = mcf10a.list, anchor.features = features)
mcf10a.combined <- IntegrateData(anchorset = mcf10a.anchors)
DefaultAssay(mcf10a.combined) <- "integrated"
all.genes <- rownames(mcf10a.combined)
mcf10a.combined <- ScaleData(mcf10a.combined, features = all.genes)
mcf10a.combined <- RunPCA(mcf10a.combined, npcs=20, features = VariableFeatures(object = mcf10a.combined))
mcf10a.combined <- RunTSNE(mcf10a.combined,dims=1:10)
FeaturePlot(mcf10a.combined,
            features = c("TFF1","KRT19","DSCAM-AS1",
                         "KRT8","KRT18","COX6C","H3F3B",
                         "PARD6B","S100A10","PSMA7","percent.mt"),reduction="tsne",slot="data")+
  FeaturePlot(mcf10a.combined,features = "dtd_FLD004",reduction="tsne",max.cutoff = 2)+
  DimPlot(mcf10a.combined,reduction="tsne",group.by = "NEP")+
  DimPlot(mcf10a.combined,reduction="tsne",group.by = "celltype")

# correlation
cancer_merge_40V <- subset(mcf10a.combined,subset=NEP=="40V" & celltype=="PC3")
exp.matrix <- t(data.frame(cancer_merge_40V[["integrated"]]@data))
var.gene <- data.frame(VariableFeatures(cancer_merge_40V))
response <- data.frame(t(data.frame(cancer_merge_40V[["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
#cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors$gene <-rownames(cors)
cors<-cors[!is.na(cors$cors),]
cors$rank <-rank(cors$cors)
pc3cor<-cors
mcf10acor<-cors
mcf7cor<-cors
mdamb231cor <- cors

gene.list <-data.frame(unique(c(mdamb231cor$gene,pc3cor$gene,mcf7cor$gene,mcf10acor$gene)))
colnames(gene.list)<-"gene"
rownames(gene.list)<-gene.list$gene
gene.list$mbamb231 <- mdamb231cor[gene.list$gene,]$cors
gene.list$pc3 <- pc3cor[gene.list$gene,]$cors
gene.list$mcf7 <- mcf7cor[gene.list$gene,]$cors
gene.list$mcf10a <- mcf10acor[gene.list$gene,]$cors

gene.list.z <- data.frame(scale(gene.list[complete.cases(gene.list),colnames(gene.list) != "gene"]))

pheatmap(gene.list[complete.cases(gene.list),colnames(gene.list) != "gene"])

gene.list.z$Mean <-rowMeans(gene.list.z)
pheatmap(gene.list.z)

ggplot(melt(gene.list.z,id.vars="Mean"),aes(x=Mean, y= value,color=variable))+geom_point()

pheatmap(gene.list.z[abs(gene.list.z$Mean)>0.0,colnames(gene.list.z) != "Mean"])
#gene.list.z[abs(gene.list.z$Mean)>1.2,]
View(cors)


#neg.cors<-subset(cors,subset=cors < -0.15)
pos.cors <- cors %>% top_n(100,cors)
neg.cors <- cors %>% top_n(-100,cors)
cors$bool <- FALSE
cors[cors$gene %in% pos.cors$gene | cors$gene %in% neg.cors$gene ,]$bool<-TRUE

mcf10a.combined <- ScaleData(mcf10a.combined, features = all.genes)
mcf10a.combined <- RunPCA(mcf10a.combined, npcs=20, features = cors[cors$bool==TRUE,]$gene)
mcf10a.combined <- RunTSNE(mcf10a.combined,dims=1:10)
FeaturePlot(subset(mcf10a.combined,subset=NEP=="40V"),
            features = c("CALM1","GNAS","KCNJ3",
                         "KRT18","KRT19","TFF1","NCOA3","nFeature_RNA"),
            reduction="tsne",
            slot="data")+
  FeaturePlot(subset(mcf10a.combined,subset=NEP=="40V"),features = "dtd_FLD004",reduction="tsne")+
  DimPlot(mcf10a.combined,reduction="tsne",group.by = "NEP")
FeatureScatter(subset(subset(mcf10a.combined,subset=NEP=="40V"),subset=celltype=="MCF10A"),feature2 = "CTSD",feature1 = "dtd_FLD004")


library(ggrepel)
ggplot(cors,aes(x=rank,y=cors))+geom_point()+geom_point(data=subset(cors,bool==TRUE),aes(x=rank,y=cors, color = "red"))+
  geom_text_repel(data=subset(cors,bool==TRUE),aes(rank,cors,label=gene),max.overlaps = 100,fontface = "italic")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")

library("org.Hs.eg.db")
library(clusterProfiler)

gene.list<-symbol2entrez(pos.cors$gene)

ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result)
View(data.frame(ego_result))

gene.plus.list<-symbol2entrez(pos.cors$gene)
gene.minus.list<-symbol2entrez(neg.cors$gene)
gene.list <- list(negative=gene.minus.list$gene_id,
                  positive=gene.plus.list$gene_id)
xx <- compareCluster(gene.list, fun="enrichKEGG")
summary(xx)
dotplot(xx)
library(pathview)

gene.list <- unlist(strsplit(xx@compareClusterResult[xx@compareClusterResult$Description=="Estrogen signaling pathway",]$geneID,"/"))
pv <- pathview(gene.data = gene.list, pathway.id = "hsa04915", species = "hsa", gene.idtype = "KEGG")
SYBOL2EG[SYBOL2EG$gene_id %in% gene.list,]

#    DimPlot(cancer_merge_40V,reduction = "tsne",group.by = "celltype")

gene.list <- unlist(strsplit(xx@compareClusterResult[xx@compareClusterResult$Description=="Proteasome",]$geneID,"/"))
pv <- pathview(gene.data = gene.list, pathway.id = "hsa03050", species = "hsa", gene.idtype = "KEGG")
SYBOL2EG[SYBOL2EG$gene_id %in% gene.list,]

mcf10a<-subset(subset(mcf10a.combined,subset=NEP=="40V"),subset=celltype=="MCF10A")
mcf10a_dtd <-data.frame(t(mcf10a[["DTD"]]@counts))
library(reshape2)
mcf10a_dtd_melt<-melt(mcf10a_dtd,measure.vars = c("FLD004","FLD010","FLD040","FLD070","FLD150","FLD500"),value.name = "Counts")
ggplot(mcf10a_dtd_melt,aes(x=variable,y=Counts))+geom_violin()+geom_jitter(size=0.1)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")

# correlation
exp.matrix <- t(data.frame(mcf10a[["RNA"]]@data))
var.gene <- data.frame(VariableFeatures(mcf10a))
response <- data.frame(t(data.frame(mcf10a[["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
#cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors$gene <-rownames(cors)
cors<-cors[!is.na(cors$cors),]
cors$rank <-rank(cors$cors)
pos.cors <- cors %>% top_n(100,cors)
neg.cors <- cors %>% top_n(-100,cors)
cors$bool <- FALSE
cors[cors$gene %in% pos.cors$gene | cors$gene %in% neg.cors$gene ,]$bool<-TRUE
ggplot(cors,aes(x=rank,y=cors))+geom_point()+geom_point(data=subset(cors,bool==TRUE),aes(x=rank,y=cors, color = "red"))+
  geom_text_repel(data=subset(cors,bool==TRUE),aes(rank,cors,label=gene),max.overlaps = 100,fontface = "italic")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+theme_bw()+
  geom_text_repel(data=subset(cors,abs(cors) > 0.25),aes(rank,cors,label=gene),max.overlaps = 100,fontface = "italic")+ 
  geom_point(data=subset(cors, abs(cors) > 0.25),aes(x=rank,y=cors, color = "red"))
FeatureScatter(subset(subset(cancer_merge1234,subset=NEP=="40V"),subset=celltype=="MCF10A"), feature1 = "CTSD", feature2 = "FLD004", group.by = "NEP")
FeatureScatter(subset(subset(cancer_merge1234,subset=NEP=="40V"),subset=celltype=="MCF10A"), feature2 = "CXCL8", feature1 = "FLD004", slot = "data")
FeatureScatter(mcf10a, feature1 = "FN1", feature2 = "FLD004", slot = "data")
DimPlot(subset(cancer_merge1234, subset=celltype=="MCF10A"), group.by = "NEP", reduction = "pca")
FeaturePlot(subset(cancer_merge1234, subset=celltype=="MCF10A"),features = "dtd_FLD004",reduction="umap")


gene.list<-symbol2entrez(pos.cors$gene)

ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result)
View(data.frame(ego_result))
gene.plus.list<-symbol2entrez(pos.cors$gene)
gene.list <- list(positive=gene.plus.list$gene_id)
xx <- compareCluster(gene.list, fun="enrichKEGG")
summary(xx)
dotplot(xx)
library(pathview)






rm(mcf10a,mcf10a_dtd,mcf10a_dtd_melt)

#
# DTD controls
#
cancer1_dtd <- Read10X(data.dir = "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX010_011/CITE-seq/AS-ele-DTD/read_count/",gene.column = 1)
cancer2_dtd <- Read10X(data.dir = "/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_10x/CITE-seq/02-ele-DTD/umi_count/",gene.column = 1)
mhsc_dtd <- Read10X(data.dir = "/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_10x/CITE-seq/03-ele-DTD/umi_count/",gene.column = 1)
