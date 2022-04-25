#https://satijalab.org/seurat/articles/integration_introduction.html
tig.list <- SplitObject(tig, split.by = "condition")
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
tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
tig.combined <- FindClusters(tig.combined, resolution = 0.1)
tig.combined <- RenameIdents(tig.combined, `0` = "TIG1-50", `1` = "TIG1-20")

DimPlot(tig, reduction = "umap",group.by = "condition")+
  DimPlot(tig, reduction = "umap")+
  DimPlot(tig.combined, reduction = "umap", group.by = "condition")+
  DimPlot(tig.combined, reduction = "umap", label = TRUE, repel = TRUE)

FeaturePlot(tig.combined,features = c("HSD17B2","ITGA8","COLEC10","DLGAP5","CENPW","dtd_FLD040"), reduction = "umap")

# remove low quality samples
tig.combined <- FindClusters(tig.combined, resolution = 0.2)
DimPlot(tig.combined)
FeatureScatter(tig.combined,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

# tig.combined <- RunPCA(tig.combined, npcs = 30, verbose = FALSE)
# tig.combined <- RunUMAP(tig.combined, reduction = "pca", dims = 1:30)
# tig.combined <- FindNeighbors(tig.combined, reduction = "pca", dims = 1:30)
# tig.combined <- FindClusters(tig.combined, resolution = 0.1)
# p1 <- DimPlot(tig.combined, reduction = "umap", group.by = "condition")
# p2 <- DimPlot(tig.combined, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2

tig.combined <- NormalizeData(tig.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
FeaturePlot(tig.combined,features=c("dtd_FLD004","TNFAIP3"))

tig.combined.nep <- subset(tig.combined,subset=condition=="NEP")
tig.combined.ctl <- subset(tig.combined,subset=condition=="CTL")


tig.marker <- FindMarkers(tig.combined.ctl,ident.1 = "TIG1-50",ident.2 = "TIG1-20",
                          thresh.use = 0.25, test.use = "wilcox", min.pct = 0.1,
                          print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf)
gene.list<-symbol2entrez(rownames(tig.marker[tig.marker$avg_log2FC>0.2 & tig.marker$p_val_adj<0.1,]))

ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE)
View(data.frame(ego_result))
aging_ego <- ego_result@result[ego_result@result$Description %in% c("aging","cell aging"),]
aging_gene<-unique(unlist(strsplit(aging_ego$geneID,"/")))
cors$bool <- FALSE

cors[cors$gene %in% gene_list_name$symbol,]$bool<-TRUE
cors[cors$gene %in% "TRPA1",]$bool<-TRUE

library(ggplot2)
library(ggrepel)
ggplot(cors,aes(x=rank,y=cors))+geom_point()+
  geom_point(data=subset(cors,bool==TRUE),aes(x=rank,y=cors, color = "red"))+
  geom_label_repel(data=subset(cors,bool==TRUE),aes(rank,cors,label=gene))
library(Seurat)
DoHeatmap(tig.bulk.seurat,features=aging_gene)

expdata<-data.frame(tig.bulk.seurat[["RNA"]]@data)
colnames(expdata)<-unlist(tig.bulk.seurat[["type_culture"]])
expdat<-expdata[rownames(expdata) %in% aging_gene,]
library(pheatmap)
pheatmap(expdat,scale="row")

library(destiny)

dm <- DiffusionMap(t(expdat))

pheatmap(expdat[,order(dm$DC1,decreasing = FALSE)],scale="row",cluster_cols = FALSE)

#tig.combined.nep <- RenameIdents(tig.combined.nep, `0` = "TIG1-50", `1` = "TIG1-20")
FeaturePlot(tig.combined.nep,features=c("dtd_FLD004","ATF3","TNFAIP3","CXCL8"))
FeaturePlot(tig.combined.ctl,features=c("dtd_FLD004","ATF3","TNFAIP3","CXCL8"))


DimPlot(tig.combined.nep)+FeaturePlot(tig.combined.nep,features=c("dtd_FLD040"))
FeaturePlot(tig.combined.nep,features=c("SAT1","FOS","FOSB","ADM"))
FeaturePlot(tig.combined.nep,features=c("TMPO","DTYMK","NCL","PCLAF"))


source("elast.CellCycle.TIG1.R")


# 
# 
# 
# tig.50 <- subset(tig.combined.nep,idents="TIG1-50")
# tig.50<- CellCycleScoring(tig.50,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
# 
# 
# tig.20 <- subset(tig.combined.nep,idents="TIG1-20")
# tig.20<- CellCycleScoring(tig.20,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
# tig.50 <- NormalizeData(tig.50,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
# tig.20 <- NormalizeData(tig.20,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
# 
# DimPlot(tig.20,reduction="umap")+FeaturePlot(tig.20,features = "dtd_FLD004")+
#   DimPlot(tig.50,reduction="umap")+FeaturePlot(tig.50,features = "dtd_FLD004")
