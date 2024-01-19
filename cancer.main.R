require(readr)
library(plyr,dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(tidyverse,R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(Seurat)
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
cancer <-cancer1
source("cancer1_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX010_011/cancer2_12/outs/filtered_feature_bc_matrix/"
cancer2 <-load.elast.data(wdir,"cancer2",100)
cancer2[["run"]]<-"second"
cancer <-cancer2
source("cancer2_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX012_013/cancer3_12/outs/filtered_feature_bc_matrix/"
cancer3 <-load.elast.data(wdir,"cancer3",100)
cancer3[["run"]]<-"third"
cancer<-cancer3
source("cancer3_Seurat.clustering.R")

wdir<- "/home/samba/sanger/shintaku/ELASTomics/20220603HiSeqX012_013/cancer4_12/outs/filtered_feature_bc_matrix/"
cancer4 <-load.elast.data(wdir,"cancer4",100)
cancer4[["run"]]<-"forth"
cancer<-cancer4
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

cancer_merge12 <- merge(cancer1,y=cancer2)
cancer_merge34 <- merge(cancer3,y=cancer4)
cancer_merge56 <- merge(cancer5,y=cancer6)
cancer_merge1234 <- merge(cancer_merge12,y=cancer_merge34)
cancer_merge123456 <- merge(cancer_merge1234,y=cancer_merge56)
rm(cancer1, cancer2, cancer3, cancer4, cancer5, cancer6, cancer_merge12, cancer_merge34, cancer_merge56, cancer_merge1234)
FeatureScatter(cancer_merge123456,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA")
cancer_merge123456 <- subset(cancer_merge123456, subset = nCount_RNA > 2000)
cancer_merge123456 <- subset(cancer_merge123456, subset = nCount_RNA < 100000)
cancer_merge123456 <- NormalizeData(cancer_merge123456, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer_merge123456  <- NormalizeData(cancer_merge123456 ,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_merge123456 <- FindVariableFeatures(cancer_merge123456, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(cancer_merge123456)
cancer_merge123456 <- ScaleData(cancer_merge123456, features = all.genes)
cancer_merge123456 <- RunPCA(cancer_merge123456, npcs=50, features = VariableFeatures(object = cancer_merge123456))
cancer_merge123456 <- JackStraw(cancer_merge123456, num.replicate = 100)
cancer_merge123456 <- ScoreJackStraw(cancer_merge123456, dims = 1:20)
JackStrawPlot(cancer_merge123456, dims = 1:20)
ElbowPlot(cancer_merge123456)
cancer_merge123456 <- RunUMAP(cancer_merge123456, dims = 1:10)
cancer_merge123456 <- RunTSNE(cancer_merge123456, npcs=50, features = VariableFeatures(object = cancer_merge123456))
DimPlot(subset(subset(cancer_merge123456, subset = NEP=="75V", invert =TRUE), subset = condition=="normal") ,reduction="umap" ,group.by = "run")
FeaturePlot(subset(subset(cancer_merge123456, subset = NEP=="75V", invert =TRUE), subset = condition=="normal") ,reduction="umap" ,features = "dtd_FLD004", max.cutoff = 4)
DimPlot(cancer_merge123456)
#
#Figure
#
VlnPlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF10A"),features = "dtd_FLD004",group.by = "condition", cols = c("#D55E00", "#0072B2"), pt.size = 0)+ scale_x_discrete(limits=c("normal", "CytoD"))+
  geom_boxplot(width = 0.1, color = "black", fill="white")  
VlnPlot(subset(subset(cancer_merge123456,subset=condition=="normal"),subset=celltype=="PC3"),features = "dtd_FLD004",group.by = "NEP", cols = c("#D55E00","#40A39A", "#0072B2"), pt.size = 0)+ scale_x_discrete(limits=c("0V", "40V", "75V"))+
  geom_boxplot(width = 0.1, color = "black", fill="white")  
RidgePlot(subset(subset(cancer_merge123456,subset=NEP=="40V"),subset=condition=="normal"), features = "dtd_FLD004",group.by = "celltype")
RidgePlot(subset(subset(subset(cancer_merge123456, subset = celltype =="MCF7"),subset=NEP=="40V"),subset=condition=="normal"), features = "dtd_FLD004",group.by = "run")
RidgePlot(subset(subset(cancer_merge123456,subset=run==c("first", "second")),subset=celltype=="MCF7"),features = "dtd_FLD004",group.by = "NEP")+
  RidgePlot(subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF10A"),features = "dtd_FLD004",group.by = "run")
FeatureScatter(subset(subset(cancer_merge123456, subset = NEP=="75V", invert =TRUE), subset = condition=="normal"),
                      feature1 = "dtd_FLD004" ,feature2 = "nCount_RNA", group.by = "celltype")
#
#PC3-voltage
#
cancer_PC3 <- subset(subset(cancer_merge123456, subset=celltype=="PC3"), subset=condition=="normal")
cancer_PC3 <- NormalizeData(cancer_PC3,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_PC3 <- NormalizeData(cancer_PC3, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer_PC3 <- ScaleData(cancer_PC3, verbose = FALSE)
cancer_PC3 <- RunPCA(cancer_PC3, npcs = 30, verbose = FALSE)
cancer_PC3 <- RunUMAP(cancer_PC3, reduction = "pca", dims = 1:30)
DimPlot(cancer_PC3,reduction="umap" ,group.by = "NEP")+
  FeaturePlot(cancer_PC3, features = "dtd_FLD004",reduction = "umap", max.cutoff = 4)
EP40V <- FindMarkers(subset(cancer_PC3, subset=NEP=="75V", invert=TRUE), group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
EP40V$lab[EP40V$p_val_adj >= 0.05] <- NA
ggplot(data=EP40V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP40V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  geom_text_repel()+theme_minimal()+theme_classic()
+xlim(-1,1)+ylim(0,50)
EP75V <- FindMarkers(subset(cancer_PC3, subset=NEP=="40V", invert=TRUE), group.by="NEP", ident.1 = "75V", min.pct = 0.25, logfc.threshold = 0)
EP75V$lab <- rownames(EP75V)
EP75V$lab[EP75V$p_val_adj >= 0.05] <- NA
ggplot(data=EP75V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP75V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  theme_minimal()+theme_classic()
rm(cancer, cancer_PC3, EP40V, EP75V)
#
# Dextran ratio
#
DTDcoc <- c(10455, 9811, 10015, 4647, 11254, 6004)
DTDmob <- c(2.26E-08,2.29E-08,2.38E-08,1.72E-08,1.44E-08,1.69E-08)
DTD40V <- subset(subset(subset(cancer_merge123456,subset=NEP=="40V"),subset=condition=="normal"), subset=celltype=="MCF10A")
DTD40V_dtd <-data.frame(t(DTD40V[["DTD"]]@counts))
DTD40V_dtd$FLD004 <- log10(DTD40V_dtd$FLD004 * (10000 / DTDcoc[1]) / (DTDmob[1] / 1E-07))
DTD40V_dtd$FLD010 <- log10(DTD40V_dtd$FLD010 * (10000 / DTDcoc[2]) / (DTDmob[2] / 1E-07))
DTD40V_dtd$FLD040 <- log10(DTD40V_dtd$FLD040 * (10000 / DTDcoc[3]) / (DTDmob[3] / 1E-07))
DTD40V_dtd$FLD070 <- log10(DTD40V_dtd$FLD070 * (10000 / DTDcoc[4]) / (DTDmob[4] / 1E-07))
DTD40V_dtd$FLD150 <- log10(DTD40V_dtd$FLD150 * (10000 / DTDcoc[5]) / (DTDmob[5] / 1E-07))
DTD40V_dtd$FLD500 <- log10(DTD40V_dtd$FLD500 * (10000 / DTDcoc[6]) / (DTDmob[6] / 1E-07))
DTD40V_dtd_melt<-melt(DTD40V_dtd,measure.vars = c("FLD004","FLD010","FLD070","FLD500"),value.name = "Counts")
ggplot(DTD40V_dtd_melt,aes(x=variable,y=Counts))+geom_violin()+geom_jitter(size=0.5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")+ylim(1,5.5)
mean(DTD40V_dtd$FLD004)
mean(DTD40V_dtd$FLD010)
mean(DTD40V_dtd$FLD070)
mean(DTD40V_dtd$FLD500)
write.table(DTD40V_dtd, file = str_c("/home/samba/public/shiomi/", "PC3", "_DTD.csv", sep = ""), append = F,sep = ",", row.names = F, quote = F)

cancer.integrated[["DTDcon"]]<-cancer.integrated[["DTD"]]
cancer.integrated[["DTDcon"]]@counts["FLD004",] <- cancer.integrated[["DTDcon"]]@counts["FLD004",]* (10000 / DTDcoc[1]) / (DTDmob[1] / 1E-07)
cancer.integrated[["DTDcon"]]@counts["FLD010",] <- cancer.integrated[["DTDcon"]]@counts["FLD010",]* (10000 / DTDcoc[2]) / (DTDmob[2] / 1E-07)
cancer.integrated[["DTDcon"]]@counts["FLD070",] <- cancer.integrated[["DTDcon"]]@counts["FLD070",]* (10000 / DTDcoc[4]) / (DTDmob[4] / 1E-07)
cancer.integrated[["DTDcon"]]@counts["FLD500",] <- cancer.integrated[["DTDcon"]]@counts["FLD500",]* (10000 / DTDcoc[6]) / (DTDmob[6] / 1E-07)
FeatureScatter(subset(subset(cancer.integrated, subset=NEP=="40V"), subset=run=="sixth", invert = TRUE), feature1 = "dtdcon_FLD004" ,feature2 = "dtdcon_FLD500", ,group.by = "celltype", slot = "count")+ylim(0,200000)

rm(DTDcoc, DTDmob, DTD40V, DTD40V_dtd, DTD40V_dtd_melt)

source("cancer.integrated.R")
cancer.integrated <- NormalizeData(cancer.integrated,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer.integrated <- NormalizeData(cancer.integrated, normalization.method = "LogNormalize", scale.factor = 1e5)
DimPlot(cancer.integrated, reduction="umap" ,group.by = "NEP")
DimPlot(subset(cancer.integrated,subset=run==c("fifth", "sixth"), invert=TRUE), reduction="umap" ,group.by = "celltype")
DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "NEP", cols = c("#D55E00", "#0072B2"))
DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "run")
DimPlot(subset(cancer.integrated,subset=NEP=="40V", invert=FALSE) ,reduction="umap" ,group.by = "run")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V") ,reduction="umap" ,features = "FLD004")
VlnPlot(subset(cancer.integrated,subset=NEP=="40V"), features = "FLD004", group.by = "celltype")
DimPlot(subset(subset(cancer.integrated,subset=run=="fifth", invert = TRUE),subset=run=="sixth", invert = TRUE) ,reduction="umap" ,group.by = "run")
VlnPlot(subset(cancer.integrated,subset=celltype=="MCF7"), features = "FLD004", group.by = "run")
VlnPlot(subset(subset(cancer.integrated,subset=NEP=="40V"),subset=run==c("fifth", "sixth"), invert = TRUE), features = "FLD004", group.by = "celltype", pt.size = 0)+
  geom_boxplot(width = 0.1, color = "black", fill="white")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "RRAD",reduction = "umap") #MCF10Atag
#DNAtag("dtd_M1AE31", "dtd_MM2E32", "dtd_P3HE33", "dtd_P3ME34", "dtd_P3NE35", "dtd_P3CE36", "dtd_MC7E37")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_M1AE31",reduction = "umap") #MCF10Atag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_MM2E32",reduction = "umap", max.cutoff = 2) #MDAMB231tag
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_P3NE35",reduction = "umap", max.cutoff = 3, min.cutoff = 1.6) #PC3tag
VlnPlot(subset(cancer.integrated,subset=run=="first"), features = "dtd_P3NE35", group.by = "celltype", slot = "data")+
  VlnPlot(subset(cancer.integrated,subset=run=="first"), features = "dtd_MC7E37", group.by = "celltype", slot = "data")
FeatureScatter(subset(cancer.integrated,subset=run=="first"), feature1 = "FLD004" ,feature2 = "dtd_MC7E37", group.by = "celltype")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), features = "dtd_MC7E37",reduction = "umap", max.cutoff = 6, min.cutoff = 2) #MCF7tag
FeatureScatter(subset(cancer.integrated,subset=celltype==c("MCF7","PC3")),feature1 = "FLD004" ,feature2 = "dtd_MC7E37")

VlnPlot(subset(cancer.integrated,subset=run==c("first", "fifth", "sixth")),features = "dtd_MC7E37",group.by = "run")
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"),
            features = c("dtd_M1AE31", "dtd_MM2E32", "dtd_P3HE33", "dtd_P3ME34", "dtd_P3NE35", "dtd_P3CE36", "dtd_MC7E37"),
            reduction = "umap", max.cutoff = 4, min.cutoff = 2)

tig <- subset(subset(cancer.integrated,subset=NEP=="0V"), subset=celltype=="PC3")
dim(tig)
summary(tig$nCount_RNA)
summary(tig$nFeature_RNA)
rm(tig)

DimPlot(subset(cancer.integrated,subset=NEP=="75V", invert=TRUE) ,reduction="umap" ,group.by = "NEP")
EP40V <- FindMarkers(subset(cancer.integrated,subset=celltype=="MCF10A", invert=TRUE), group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
EP40V$lab[EP40V$p_val_adj >= 0.05] <- NA
ggplot(data=EP40V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP40V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  geom_text_repel()+theme_minimal()+theme_classic()
ggplot(EP40V, aes(x=avg_log2FC)) + geom_histogram(color="black", fill="white")+scale_y_log10()

FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "umap", features = "CYBA") #MCF7
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "umap", features = "S100A2") #PC-3
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "umap", features = "CAVIN3") #MDAMB231
FeaturePlot(subset(cancer.integrated,subset=NEP=="40V"), reduction = "umap", features = "DUSP1") #MCF10A

source("cancer.MCF10A.R")
#Tukey's test
FACS <- subset(subset(cancer_merge123456,subset=run==c("fifth", "sixth")),subset=celltype=="MCF10A")
FACS1 <-cbind(data.frame(FACS[["DTD"]]@counts["FLD004",]), data.frame(FACS[["condition"]]))
colnames(FACS1) <- c("FLD004", "Hash")
t.test(FLD004 ~ Hash, data = FACS1)

FACS <- subset(subset(cancer_merge123456,subset=condition=="normal"),subset=celltype=="PC3")
FACS1 <-cbind(data.frame(FACS[["DTD"]]@counts["FLD004",]), data.frame(FACS[["NEP"]]))
colnames(FACS1) <- c("FLD004", "Cell")
amod <- aov(FLD004 ~ Cell, data = FACS1)
TukeyHSD(amod)

FACS <- subset(cancer.integrated,subset=NEP=="40V")
FACS1 <-cbind(data.frame(FACS[["DTD"]]@counts["FLD004",]), data.frame(FACS[["celltype"]]))
colnames(FACS1) <- c("FLD004", "Cell")
amod <- aov(FLD004 ~ Cell, data = FACS1)
TukeyHSD(amod)

rm(FACS, FACS1)
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
cancer.integrated<-CellCycleScoring(cancer.integrated,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
VlnPlot(subset(subset(cancer.integrated,subset=NEP=="40V"), subset=condition=="normal"),features = "dtd_FLD004")
RidgePlot(subset(cancer.integrated,subset=NEP=="40V"),features = "dtd_FLD004")
FeatureScatter(subset(cancer.integrated,subset=NEP=="40V"),feature2 = "dtd_FLD004",feature1 = "S.Score")
FeatureScatter(subset(cancer.integrated,subset=NEP=="40V"),feature2 = "dtd_FLD004",feature1 = "G2M.Score")
FeatureScatter(subset(cancer.integrated,subset=NEP=="40V"),feature1 = "G2M.Score" ,feature2 = "S.Score")
FeatureScatter(subset(subset(cancer.integrated,subset=NEP=="40V"), subset=celltype=="MCF7"),feature2 = "dtd_FLD004",feature1 = "S.Score")
FeatureScatter(subset(subset(cancer.integrated,subset=NEP=="40V"), subset=celltype=="MCF7"),feature2 = "dtd_FLD004",feature1 = "G2M.Score")
#
#mcf7-MDAMB231 subset 
#
MCFMDA <- subset(subset(subset(cancer.integrated, subset=celltype=="PC3", invert=TRUE), subset=celltype=="MCF10A", invert=TRUE), subset=NEP=="40V")
FeatureScatter(MCFMDA,feature1 = "nCount_RNA" ,feature2 = "FLD004", group.by = "run")
VlnPlot(subset(MCFMDA, subset =celltype=="MCF7"),feature = "nCount_RNA", group.by = "run")
FeatureScatter(MCFMDA,feature1 = "nCount_RNA" ,feature2 = "nFeature_RNA", group.by = "celltype")
FeatureScatter(MCFMDA,feature1 = "G2M.Score" ,feature2 = "S.Score")
Idents(MCFMDA)<-MCFMDA[["NEP"]]
MCFMDA<-NormalizeData(MCFMDA,assay="RNA",normalization.method = "LogNormalize", scale.factor = 1e5)
MCFMDA<-NormalizeData(MCFMDA,assay="DTD",normalization.method = "CLR", scale.factor = 1e2)
#
#Correlation
#
library(glmnet)
exp.matrix <- t(data.frame(MCFMDA[["RNA"]]@data))
var.gene <- data.frame(VariableFeatures(MCFMDA))
response <- data.frame(t(data.frame(MCFMDA[["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
cors.sig<-subset(subset(cors,subset=pval<1e-4),subset=abs(cors) > 0.2)
cors.pos<-subset(subset(cors,subset=pval<1e-4),subset=cors > 0.2)
cors.neg<-subset(subset(cors,subset=pval<1e-4),subset=cors < -0.2)
ego_gene <- c("BIRC3", "CXCL1", "CXCL2", "CXCL3", "IL6", "IRF1", "JUNB", "LIF", "NFKBIA", "CCL20", "TNFAIP3")
cors$bool <- FALSE
cors[cors$gene %in% ego_gene,]$bool<-TRUE
ggplot(data=subset(cors, abs(cors) <= 0.15),aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.15),aes(x=rank,y=cors, color = "red"))+NoLegend()+theme_classic()+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1)
cors <- cors[order(cors$rank, decreasing=F),]

library(gplots)
MCF7 <- FindMarkers(subset(subset(subset(cancer.integrated, subset=celltype=="PC3", invert=TRUE), subset=celltype=="MCF10A", invert=TRUE), subset=NEP=="0V"),
                           group.by="celltype", ident.1 = "MCF7", ident.2 = "MDAMB231", min.pct = 0.25, logfc.threshold = 0)
MCF7$gene <- rownames(MCF7)
MCF7 <- inner_join(MCF7, cors, by="gene")
rownames(MCF7) <- MCF7$gene 
ggplot(MCF7,aes(x=avg_log2FC, y=cors))+geom_point(color = "black")+theme_classic()
ggplot(MCF7,aes(x=avg_log2FC,y=cors))+geom_point()+
  geom_point(data=subset(MCF7, abs(avg_log2FC) > 2.5),aes(x=avg_log2FC,y=cors, color = "red"))+xlim(-10, 10)+theme_classic()

data <- list(MDA = rownames(MCF7[MCF7$avg_log2FC < -2, ]), 
             MCF7 = rownames(MCF7[MCF7$avg_log2FC > 2, ]), 
             posi = cors.pos$gene, nega = cors.neg$gene)
venn(data)
source("cancer.GSEA.R")
#
#egoGo
#
library(clusterProfiler)
library(org.Hs.eg.db)
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}
gene <- symbol2entrez(rownames(cors.neg))
gene <- symbol2entrez(intersect(rownames(cors.neg), rownames(MCF7[MCF7$avg_log2FC < -2, ])))
gene <- symbol2entrez(setdiff(rownames(cors.neg), rownames(MCF7[MCF7$avg_log2FC < -2, ])))
gene <- symbol2entrez(intersect(rownames(cors.pos), rownames(MCF7[MCF7$avg_log2FC > 2, ])))
gene <- symbol2entrez(setdiff(rownames(cors.pos), rownames(MCF7[MCF7$avg_log2FC > 2, ])))
gene <- symbol2entrez(rownames(MCF7[MCF7$cors < -0.1 & MCF7$avg_log2FC > -2.5 & MCF7$avg_log2FC < 2.5,]))
gene <- symbol2entrez(rownames(MCF7[MCF7$cors < -0.1 & MCF7$avg_log2FC < -2.5,]))
gene <- symbol2entrez(rownames(MCF7[MCF7$cors < -0.15,]))
gene <- symbol2entrez(rownames(MCF7[MCF7$cors > 0.2,]))
gene <- symbol2entrez(rownames(MCF7[MCF7$cors > 0.2 & MCF7$avg_log2FC > -2.5,]))
ego_result <- enrichGO(gene          = gene$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
barplot(ego_result, drop=TRUE, showCategory=30)+theme_classic()
cnetplot(ego_result)
View(data.frame(ego_result))
ego_gene <- as.vector(str_split(as.data.frame(ego_result)["GO:0045296", "geneID"], pattern = "/")[[1]])

kegg_result <- enrichKEGG(gene=gene$gene_id,
                          pvalueCutoff=0.05)
barplot(kegg_result, drop=TRUE, showCategory=30)
kegg_result <- pairwise_termsim(kegg_result)
emapplot(kegg_result)


# elastic net_mcf10a
exp.matrix <- t(data.frame(mcf10a[["RNA"]]@data))
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
alpha <- seq(0, 1, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = best.alpha)
best.lambda <- m$lambda.min
en.model <- glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian",
                   lambda = best.lambda, alpha = best.alpha)
en.model.beta <-data.frame(coef(en.model, s="lambda.min"))
en.model.beta$gene <- rownames(en.model.beta)
en.model.beta$cors <- cors.sig[en.model.beta$gene,]$cors
en.model.beta <- na.omit(en.model.beta)
ggplot(en.model.beta,aes(x=cors,y=s1))+geom_point()+
  geom_point(data=subset(en.model.beta,gene %in% rownames(en.model.beta[abs(en.model.beta$s1)>0.1,])),aes(y=s1,x=cors,color="red"))+
  geom_text_repel(data=subset(en.model.beta,gene %in% rownames(en.model.beta[abs(en.model.beta$s1)>0.1,])),aes(y=s1,x=cors,label=gene,fontface = "italic"))+
  theme_bw()+NoLegend()+theme_bw()+NoLegend()

FeatureScatter(mcf10a,feature1 = "FLD004" ,feature2 = "HSP90AB1", slot = "data", group.by = "celltype", cols = c("#79A72C", "#1BB4B8"))
FeatureScatter(mcf10a,feature1 = "FLD004" ,feature2 = "DSCAM-AS1", slot = "data", group.by = "celltype")
VlnPlot(mcf10a, features = "HSP90AB1", group.by = "celltype")

