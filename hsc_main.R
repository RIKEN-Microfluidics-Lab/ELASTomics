require(readr)
library(plyr)
library(dplyr,tidyr)
library(ggplot2)
library(tidyverse,R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(Seurat)
library(seqinr)
library(stringr)
rdir <- "/home/samba/public/shintaku/github/ELASTomics/"

# load ELASTomics data from an output of cellranger with cite-seq pipeline
source("elast.load.elast.data.R")

# control data
wdir <- "/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_10x/mhsc_1_14c_2/outs/filtered_feature_bc_matrix/"
tig <-load.elast.data(wdir,"NEP-",100)
tig <- NormalizeData(tig, normalization.method = "LogNormalize", scale.factor = 1e5)
# mitochondrial gene percent and remove dead cells
tig[["percent.mt"]] <- PercentageFeatureSet(tig, pattern = "^mt-")
tig <- subset(tig, subset= percent.mt<5)
tig <- subset(tig,subset=nCount_RNA >1000)

tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 300)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tig), 40)
# plot variable features with labels
plot1 <- VariableFeaturePlot(tig)
plot1 <- LabelPoints(plot = tig, points = top10)
plot1

# normalize the dtd data with "RC" option.
tig <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig <- NormalizeData(tig,assay="ADT",normalization.method = "CLR",scale.factor = 1e2)
#
# PCA and visualize
all.genes <- rownames(tig)
tig <- ScaleData(tig, features = all.genes)
tig <- RunPCA(tig, npcs=20, features = VariableFeatures(object = tig))
DimPlot(tig, reduction = "pca")
source("Seurat.clustering.R")
JackStrawPlot(tig, dims = 1:20)
ElbowPlot(tig)

#Celltype
tig <- FindNeighbors(tig, dims = 1:18)
tig <- FindClusters(tig, resolution = 0.8)
DimPlot(tig, label = TRUE)
new.cluster.ids <- c("Gra", "Gra", "Gra","Gra", "3GMP", "Mono","Gra", "3EryP/MEP", "1MPP","4Erythroid", "2CMP", "Bcell","Bcell", "5Erythrocyte", "CLP","Gra", "2CMP", "Other","Other")
names(new.cluster.ids) <- levels(tig)
tig <- RenameIdents(tig, new.cluster.ids)
tig[["Celltype"]] <- Idents(tig)

#CMP or MkP
tig1 <- subset(tig, subset=Celltype=="2CMP")
tig1 <- FindNeighbors(tig1, dims = 1:10)
tig1 <- FindClusters(tig1, resolution = 0.2)
DimPlot(tig1, label = FALSE)
new.cluster.ids <- c("MkP", "2CMP", "2CMP","BaP")
names(new.cluster.ids) <- levels(tig1)
tig1 <- RenameIdents(tig1, new.cluster.ids)
tig1[["Celltype"]] <- Idents(tig1)

#HSC or MPP
tig2 <- subset(tig, subset=Celltype=="1MPP")
tig2 <- FindNeighbors(tig2, dims = 1:20)
tig2 <- FindClusters(tig2, resolution = 0.8)
DimPlot(tig2, label = FALSE)
new.cluster.ids <- c("0HSC", "1MPP", "1MPP", "1MPP", "1MPP")
names(new.cluster.ids) <- levels(tig2)
tig2 <- RenameIdents(tig2, new.cluster.ids)
tig2[["Celltype"]] <- Idents(tig2)

#Reconstruction
tig3 <- subset(tig, idents = c("2CMP", "1MPP"), invert = TRUE)
tig12 <- merge(tig1, y=tig2)
tig4 <- merge(tig3, y=tig12)
tig4 <- FindVariableFeatures(tig4, selection.method = "vst", nfeatures = 300)
all.genes <- rownames(tig4)
tig4 <- ScaleData(tig4, features = all.genes)
tig4 <- RunPCA(tig4, npcs=20, features = VariableFeatures(object = tig))
tig4 <- RunUMAP(tig4, dims = 1:19)
DimPlot(tig4, reduction = "umap", label =TRUE, pt.size = 0.5, group.by = "Celltype") + NoLegend()
tig <- tig4
tig <- subset(tig,subset=nCount_RNA < 60000)
#rm(tig1, tig2, tig12, tig3, tig4)
dim(tig)
summary(tig$nCount_RNA)
summary(tig$nFeature_RNA)
as.data.frame(tig[["Celltype"]]) %>% dplyr::count(Celltype) 

#Gene expression
FeaturePlot(tig,features = c("Camp", "Chil3","Mmp9"),reduction="umap", ncol = 3) #Granulocyte
FeaturePlot(tig,features = c("Cd48"),reduction="umap") #Granulocyte
FeaturePlot(tig,features = c("Elane", "Vcam1","Mpo"),reduction="umap", ncol = 3) #GMP
FeaturePlot(tig,features = c("Ccr2", "Klf4","Irf8"),reduction="umap", ncol = 3) #Mono
FeaturePlot(tig,features = c("Cox6a2", "Dntt","Il7r"),reduction="umap", ncol = 3) #CLP
FeaturePlot(tig,features = c("Ly6a", "Hlf", "Rgs1"),reduction="umap", ncol = 3) #MPP+HSC
VlnPlot(tig, features = c("Hlf"), log = TRUE,  slot = "count") #HSC
FeaturePlot(tig,features = c("Cd79a", "Vpreb1", "Pax5"),reduction="umap", ncol = 3) #other(B-cell)
FeaturePlot(tig,features = c("Prss34"),reduction="umap") #BaP
FeaturePlot(tig,features = c("Itga2b", "Gata2", "Pf4"),reduction="umap", ncol = 3) #CMP+MkP
VlnPlot(tig, features = c("Pf4"), log = TRUE,  slot = "count") #MkP
FeaturePlot(tig,features = c("Klf1", "Mt2", "Gata1"),reduction="umap", ncol = 3) #EryP/CMP
FeaturePlot(tig,features = c("Gypa", "Slc4a1", "Hba-a1"),reduction="umap", ncol = 3) #Erythroid
FeaturePlot(tig,features = "percent.mt",reduction="umap")
Fgene <- c("Cd79a", "Vpreb1", "Pax5", "Cox6a2", "Dntt","Il7r", "Itga2b", "Gata2", "Pf4", "Klf1", "Mt2", "Gata1", "Gypa", "Slc4a1", "Hba-a1", "Elane", "Vcam1","Mpo", "Camp", "Chil3","Mmp9", "Ly6a", "Hlf", "Rgs1", "Ccr2", "Klf4","Irf8")
DoHeatmap(tig, group.by = "Celltype", features = Fgene)


#CITE-seq
tig <- subset(tig, subset= adt_ADT130<6)
tig <- subset(tig, subset= adt_ADT429<6)
FeaturePlot(tig,features="adt_ADT130",max.cutoff = 3) #ADT130 = Ly6a(Sca-1)
VlnPlot(tig, features = c("adt_ADT130"), log = TRUE,  slot = "count") 
FeaturePlot(tig,features="adt_ADT429",max.cutoff = 2.5) #ADT429 = CD48
FeaturePlot(tig,features="Cd48") #ADT429 = CD48
VlnPlot(tig, features = c("adt_ADT429"), log = TRUE,  slot = "count")
FeaturePlot(tig,features="adt_ADT012",max.cutoff = 200) #ADT012 = c-kit(CD117)
VlnPlot(tig, features = c("adt_ADT012"), log = TRUE,  slot = "count")
FeaturePlot(tig,features="adt_ADT203",max.cutoff = 2.5) #ADT203 = CD150(SLAM)
VlnPlot(tig, features = c("adt_ADT203"), log = TRUE,  slot = "count")

#ELASTomics
FeaturePlot(tig,features="FLD004",max.cutoff = 15) 
VlnPlot(tig, features = c("dtd_FLD004"), log = TRUE,  slot = "data")
p <- pivot_longer(as.data.frame(t(tig[["DTD"]]@data)), cols = 1:6, names_to = "Categories", values_to = "Values")
g <- ggplot(p, aes(x = Categories, y = Values))
g+geom_violin()+scale_y_log10()+geom_jitter(width=0.3)
HSC_dtd <- cbind(as.data.frame(tig[[c("Celltype")]]), as.data.frame(t(tig[["DTD"]]@data)))
TukeyHSD(aov(FLD500~Celltype, data = HSC_dtd))

#ADT-DTD
tig1 <- NormalizeData(tig,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig1 <- NormalizeData(tig1,assay="ADT",normalization.method = "CLR",scale.factor = 1e2)
tig1 <- subset(tig1, idents = c("Other", "Bcell"), invert= TRUE)
#HSC_dtd <- cbind(as.data.frame(tig1[[c("Celltype")]]), as.data.frame(t(tig1[["DTD"]]@data)))
#HSC_dtd <- 100 * (table(HSC_dtd[HSC_dtd$FLD004 > 1, ]$Celltype) / table(HSC_dtd$Celltype))
tig1 <- subset(tig1, subset= dtd_FLD004>1)
VlnPlot(object = tig1, features='dtd_FLD004')
#HSC_dtd <- cbind(as.data.frame(tig1[[c("Celltype")]]), as.data.frame(t(tig1[["DTD"]]@data)))
tig1 <- subset(tig1, idents = c("0HSC", "1MPP", "2CMP", "3EryP/MEP", "4Erythroid", "5Erythrocyte"))  #Erythroid
#tig1 <- subset(tig1, idents = c("0HSC", "1MPP","2CMP",  "3GMP", "Gra", "Mono", "BaP"))               #Granulocyte

FeatureScatter(object = tig1, feature1 = 'adt_ADT130', feature2 = 'dtd_FLD500',  slot = "count")
FeatureScatter(object = tig1, feature1 = 'Sptb', feature2 = 'dtd_FLD004')
VlnPlot(tig1, features = c("adt_ADT429"), slot = "data")
VlnPlot(tig1, features = c("Sptb"), slot = "data")
VlnPlot(tig1, features = c("percent.mt"), slot = "data")
VlnPlot(tig1, features = "dtd_FLD004", slot = "data")
FeaturePlot(tig1,features="dtd_FLD500",max.cutoff = 5) 
tig1 <- subset(tig1, subset= dtd_FLD004>1)
VlnPlot(tig1, features = "dtd_FLD004", pt.size = 0.5, slot = "data")
FeatureScatter(tig1, feature2 = "Spta1",feature1 = "FLD004")

#Correlation
library(glmnet)
exp.matrix <- t(data.frame(tig1[["RNA"]]@data))
var.gene <- data.frame(VariableFeatures(tig1))
response <- data.frame(t(data.frame(tig1[["DTD"]]@data)))
rownames(response) <- rownames(exp.matrix)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
cors$bool <- FALSE
cors[cors$gene %in% rownames(cors %>% top_n(5, cors)),]$bool<-TRUE
cors[cors$gene %in% rownames(cors %>% top_n(-5, cors)),]$bool<-TRUE
cors.sig<-subset(subset(cors,subset=pval<1e-4),subset=abs(cors)>0.2)
cors.pos<-subset(subset(cors,subset=pval<1e-4),subset = cors > 0.2)
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.2),aes(x=rank,y=cors, color = "red"))+theme_classic() +ylim(-0.35, 0.35) +NoLegend()+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1) 

#egoGo
library(clusterProfiler)
library(org.Mm.eg.db)
#symbol2entrez<-function(gene.symbol){
#  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
#  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
#  gene.list<-gene.list[!is.na(gene.list$gene_id),]
#  return(gene.list)
#}
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[gene.list$symbol,1] #correlation
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
gene_list_log2fc_cor <- symbol2entrez.order(cors)
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Mm.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
View(gse_result_cor@result)
ridgeplot(gse_result_cor,showCategory = 70)

g1 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0030507") #spectrin binding
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0016887") #ATP hydrolysis activity
g3 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0003779") #actin binding
g4 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0140657") #ATP-dependent activity
gridExtra::grid.arrange(g3, g4, nrow = 2) 




MCF7 <- FindMarkers(tig1, group.by="Celltype", ident.1 = c("0HSC", "1MPP", "2CMP"), ident.2 = c("3EryP/MEP", "4Erythroid", "5Erythrocyte"), min.pct = 0.25, logfc.threshold = 0)
MCF7$gene <- rownames(MCF7)
rownames(MCF7) <- MCF7$gene 
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[gene.list$symbol,2] #foldchange
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
gene_list_log2fc_fol <- symbol2entrez.order(MCF7)  
gse_result_fol<- gseGO(geneList     = gene_list_log2fc_fol,
                       OrgDb        = org.Mm.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
View(gse_result_fol@result)


dotplot(gse_result_cor,showCategory = 12)
dotplot(gse_result_fol,showCategory = 12)
GSEA_cor <- gse_result_cor@result
GSEA_cor <- GSEA_cor[GSEA_cor$p.adjust < 0.05, ]
GSEA_cor <- GSEA_cor[,c(1, 2, 4)]
colnames(GSEA_cor) <- c("ID", "Description", "Cor_Score")
GSEA_fol <- gse_result_fol@result
GSEA_fol <- GSEA_fol[GSEA_fol$p.adjust < 0.05, ]
GSEA_fol <- GSEA_fol[,c(1, 2, 4)]
colnames(GSEA_fol) <- c("ID", "Description", "Fol_Score")
GSEA_full <- full_join(x = GSEA_cor, y= GSEA_fol, by = "ID")
GSEA_full[is.na(GSEA_full)] <- 0
GSEA_full[GSEA_full$Description.x == 0,2] <- GSEA_full[GSEA_full$Description.x == 0,4]
rownames(GSEA_full) <- GSEA_full$ID
GSEA_full <- GSEA_full[,c(2, 3, 5)]
for (d in 1:nrow(GSEA_full)){
  name <- rownames(GSEA_full)[d]
  GSEA_full[d,] <- c(GSEA_full[d,1], gse_result_cor@result[name, 4], gse_result_fol@result[name, 4])
}
rm(d, name)
GSEA_full[is.na(GSEA_full)] <- 0
rownames(GSEA_full) <- GSEA_full$Description.x
GSEA_full <- GSEA_full[,c(2, 3)]
GSEA_full[, 1] <- as.numeric(GSEA_full[, 1])
GSEA_full[, 2] <- as.numeric(GSEA_full[, 2])
GSEA_full <- as.matrix(GSEA_full)
pheatmap::pheatmap(t(GSEA_full))
#
#enrichGo
#
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}
gene <- symbol2entrez(rownames(cors.sig))
ego_result <- enrichGO(gene          = gene$gene_id, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
barplot(ego_result, drop=TRUE)




#
ego_gene <- c("Ank1", "Cpox", "Spta1", "Sptb", "Uros", "Blvrb", "Actg1", "Arhgdib", "Arpc1b", "Tmsb10", "S100a10", "Gmfg", "Arpc2")
cors$bool <- FALSE
cors[cors$gene %in% ego_gene,]$bool<-TRUE
ggplot(data=subset(cors, abs(cors) <= 0.25),aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.25),aes(x=rank,y=cors, color = "red"))+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1) +theme_bw() +NoLegend()




# heatmap of correlated genes
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.1)
exp.matrix <- data.frame(tig1[["RNA"]]@data)
dtd <-FetchData(object = tig1, vars = c("dtd_FLD500"))
rownames(dtd)<-colnames(exp.matrix)
out<-pheatmap(exp.matrix[rownames(cors.sig),order(dtd$dtd_FLD500,decreasing = FALSE)],
              cluster_cols = FALSE, show_colnames = FALSE,
              annotation_col = dtd,
              clustering_distance_rows = "correlation")

# elastic net
exp.matrix <- t(data.frame(tig1[["RNA"]]@data))
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
ggplot(en.model.beta, aes(x = cors, y = s1, label = gene))+geom_point()+theme_bw() + #geom_text_repel()
  geom_text_repel(aes(label=ifelse(abs(s1) > 0.1,as.character(gene),'')))




#
# re-scale the dtd data with the concentration of the dtd molecules in the solution
# 
source("elast.rescale.dtd.R")
# FLD004,FLD010,FLD040,FLD070,FLD150,FLD500, and others
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "RC",scale.factor = 1e2)
concentration<- data.frame(c(2,6,9,3,9,3,1,1,1,1,1,1,1))
tig.combined.nep<-rescale.dtd.data(tig.combined.nep,concentration)

tig.dtd.scale<-tig.combined.nep[["DTD"]]@data
#tig.dtd.scale <- tig.dtd.scale[c("FLD004","FLD010","FLD070","FLD500"),]


#
# compute a radius for a single case
source("elast.comp.radii.R")
#
# S.radii: input Stokes radii of DTD
# tig.dtd.scale; scaled amount of DTD imported to cells
# The nrow of the tig.dtd.scale must match with the length of S.radii
#
S.radii <- c(1.4e-9, 2.3e-9, 4.5e-9, 6.3e-9, 8.5e-9, 15.1e-9) # Stokes radii of DTD

#S.radii <- c(1.4e-9, 2.4e-9, 1.4e-9, 1.4e-9, 1.4e-9, 1.4e-9)
#fm <-elast.comp.radius(tig.dtd.scale[1:4,"NEP-CCCTTAGGTCAAACGG"],S.radii,TRUE)
#
# compute radii for multiple cases
res<-elast.comp.radii(tig.dtd.scale[1:4,],S.radii,TRUE)
Res <- res[[1]]
dtd.predict<-res[[2]]
#
# integrate elast result into Seurat object
# 

source("elast.integrated.R")

#
# extract explanatory variables via glmnet
#
source("elast.glmnet.R")
en.model.nonzero.beta <- subset(en.model.beta, subset=abs(s0)>0.001)
#source("elast.biomaRt.R")
en.model.plus.beta <- subset(en.model.beta, subset=s0> 0.001)
en.model.minus.beta <- subset(en.model.beta, subset=s0< -0.001)

source(file.path(rdir,"elast.integrated.cor.visualize.R"))

source(file.path(rdir,"elast.integrated.cor.R"))

source(file.path(rdir,"elast.integrated.cor.visualize.R"))

source(file.path(rdir,"elast.integrated.cor.diffusionmap.R"))

 source("elast.integrate.data.R")


