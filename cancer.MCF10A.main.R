library(ggplot2)
library(ggrepel)
library(Seurat)
library(pheatmap)
library(destiny)

tig <- subset(cancer,subset=celltype=="MCF10A")
RidgePlot(tig, features = "dtd_FLD004",group.by = "NEP") 

#https://satijalab.org/seurat/articles/integration_introduction.html
tig.list <- SplitObject(tig, split.by = "NEP")
tig.list <- lapply(X = tig.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = tig.list)

tig.anchors <- FindIntegrationAnchors(object.list = tig.list, anchor.features = features)
tig.combined <- IntegrateData(anchorset = tig.anchors)
DefaultAssay(tig.combined) <- "integrated"
tig.combined <- ScaleData(tig.combined, verbose = FALSE)
tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)

DimPlot(tig.combined, reduction = "pca",group.by = "NEP")
FeaturePlot(tig.combined,features = "dtd_FLD004", reduction = "pca")
tig.combined <- NormalizeData(tig.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
tig1 <- subset(tig.combined,subset=NEP=="40V")
tig1 <- subset(tig,subset=NEP=="40V")

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
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.25)
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+ylim(-0.3, 0.3) +theme_bw()
geom_point(data=subset(cors, abs(cors) > 0.25),aes(x=rank,y=cors, color = "red"))
FeatureScatter(tig1, feature2 = "dtd_FLD010",feature1 = "dtd_FLD500")

#egoGo
library(clusterProfiler)
library(org.Mm.eg.db)
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Mm.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}
gene <- symbol2entrez(append(rownames(cors[cors$cors > 0.25,]),rownames(cors[cors$cors < -0.25,])))
#gene <- symbol2entrez(rownames(cors.sig))
ego_result <- enrichGO(gene          = gene$gene_id, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))
ego_gene <- c("Ank1", "Cpox", "Spta1", "Sptb", "Uros", "Blvrb", "Actg1", "Arhgdib", "Arpc1b", "Tmsb10", "S100a10", "Gmfg", "Arpc2")
cors$bool <- FALSE
cors[cors$gene %in% ego_gene,]$bool<-TRUE
ggplot(data=subset(cors, abs(cors) <= 0.25),aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.25),aes(x=rank,y=cors, color = "red"))+
  geom_text_repel(data=subset(cors,bool==TRUE), aes(rank,cors,label=gene), min.segment.length = 0.1) +theme_bw() +NoLegend()
