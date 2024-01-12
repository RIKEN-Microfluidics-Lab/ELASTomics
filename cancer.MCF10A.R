#
#MCF10A lot
#
cancer_MCF10 <- subset(subset(cancer_merge123456, subset=celltype=="MCF10A"), subset=condition=="normal")
cancer_MCF10 <- subset(subset(cancer.integrated, subset=celltype=="MCF10A"), subset=run=="fifth", invert=TRUE)
cancer_MCF10 <- NormalizeData(cancer_MCF10,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_MCF10 <- NormalizeData(cancer_MCF10, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer_MCF10 <- RunUMAP(cancer_MCF10, dims = 1:10)
DimPlot(cancer_MCF10, reduction="umap" ,group.by = "run")
DimPlot(cancer_MCF10, reduction="umap" ,group.by = "NEP")
VlnPlot(cancer_MCF10, features = "FLD004", group.by = "run", slot = "data")
FeaturePlot(cancer_MCF10, features = "FLD004",reduction = "umap")
cancer.list <- SplitObject(cancer_MCF10, split.by = "run")
cancer.list <- lapply(X = cancer.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = cancer.list)
cancer.anchors <- FindIntegrationAnchors(object.list = cancer.list, anchor.features = features)
cancer_MCF10.combined <- IntegrateData(anchorset = cancer.anchors)
DefaultAssay(cancer_MCF10.combined) <- "integrated"
cancer_MCF10.combined <- ScaleData(cancer_MCF10.combined, verbose = FALSE)
cancer_MCF10.combined <- NormalizeData(cancer_MCF10.combined,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_MCF10.combined <- NormalizeData(cancer_MCF10.combined, normalization.method = "LogNormalize", scale.factor = 1e5)
cancer_MCF10.combined <- RunPCA(cancer_MCF10.combined, npcs = 30, verbose = FALSE)
cancer_MCF10.combined <- RunUMAP(cancer_MCF10.combined, reduction = "pca", dims = 1:30)
DimPlot(cancer_MCF10.combined,reduction="umap" ,group.by = "NEP", cols = c("#D55E00", "#0072B2"))+
  DimPlot(cancer_MCF10.combined,reduction="umap" ,group.by = "run")+
  FeaturePlot(cancer_MCF10.combined, features = "dtd_FLD004",reduction = "umap")
VlnPlot(cancer_MCF10.combined, features = "FLD004", group.by = "run")

EP40V <- FindMarkers(cancer_MCF10, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
EP40V$lab[EP40V$p_val_adj >= 0.05] <- NA
ggplot(data=EP40V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP40V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  geom_text_repel()+theme_minimal()+theme_classic()
#
#GSEA
#
library("org.Hs.eg.db")
library(clusterProfiler)
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[gene.list$symbol,2] #LogFoldchange
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
gene_list_log2fc_cor <- symbol2entrez.order(EP40V)  
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "ALL",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
View(gse_result_cor@result)
dotplot(gse_result_cor,showCategory = 12)
ridgeplot(gse_result_cor,showCategory = 12)

g1 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0008219") #cell death
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0033554") #cellular response to stress
gridExtra::grid.arrange(g1, g2, nrow = 2) 


GSEA_MCF10A <- gse_result_cor@result
GSEA_MCF10A <- GSEA_MCF10A[GSEA_MCF10A$p.adjust < 0.05, ]
GSEA_MCF10A <- GSEA_MCF10A[,c(1, 2, 4)]
colnames(GSEA_MCF10A) <- c("ID", "Description", "Cor_Score")
rm(cancer_MCF10)

cancer_MCF7 <- subset(cancer.integrated, subset=celltype=="MCF7")
cancer_MCF7 <- NormalizeData(cancer_MCF7,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_MCF7 <- NormalizeData(cancer_MCF7, normalization.method = "LogNormalize", scale.factor = 1e5)
EP40V <- FindMarkers(cancer_MCF7, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
gene_list_log2fc_cor <- symbol2entrez.order(EP40V)  
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       verbose      = FALSE)
View(gse_result_cor@result)
dotplot(gse_result_cor,showCategory = 12)
GSEA_MCF7 <- gse_result_cor@result
GSEA_MCF7 <- GSEA_MCF7[GSEA_MCF7$p.adjust < 0.05, ]
GSEA_MCF7 <- GSEA_MCF7[,c(1, 2, 4)]
colnames(GSEA_MCF7) <- c("ID", "Description", "Cor_Score")
rm(cancer_MCF7)

cancer_PC3 <- subset(cancer.integrated, subset=celltype=="PC3")
cancer_PC3 <- NormalizeData(cancer_PC3,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_PC3 <- NormalizeData(cancer_PC3, normalization.method = "LogNormalize", scale.factor = 1e5)
EP40V <- FindMarkers(cancer_PC3, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
gene_list_log2fc_cor <- symbol2entrez.order(EP40V)  
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       verbose      = FALSE)
View(gse_result_cor@result)
dotplot(gse_result_cor,showCategory = 12)
GSEA_PC3 <- gse_result_cor@result
GSEA_PC3 <- GSEA_PC3[GSEA_PC3$p.adjust < 0.05, ]
GSEA_PC3 <- GSEA_PC3[,c(1, 2, 4)]
colnames(GSEA_PC3) <- c("ID", "Description", "Cor_Score")
rm(cancer_PC3)

cancer_MDA <- subset(cancer.integrated, subset=celltype=="MDAMB231")
cancer_MDA <- NormalizeData(cancer_MDA,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_MDA <- NormalizeData(cancer_MDA, normalization.method = "LogNormalize", scale.factor = 1e5)
EP40V <- FindMarkers(cancer_MDA, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V$lab <- rownames(EP40V)
gene_list_log2fc_cor <- symbol2entrez.order(EP40V)  
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "MF",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       verbose      = FALSE)
View(gse_result_cor@result)
dotplot(gse_result_cor,showCategory = 12)
GSEA_MDA <- gse_result_cor@result
GSEA_MDA <- GSEA_MDA[GSEA_MDA$p.adjust < 0.05, ]
GSEA_MDA <- GSEA_MDA[,c(1, 2, 4)]
colnames(GSEA_MDA) <- c("ID", "Description", "Cor_Score")
rm(cancer_MDA)
rm(EP40V, gene_list_log2fc_cor, gse_result_cor)
#
#GSEA
#
cancer_MCF10 <- subset(subset(cancer.integrated, subset=celltype=="MCF10A"), subset=run=="fifth", invert=TRUE)
cancer_MCF10 <- NormalizeData(cancer_MCF10,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_MCF10 <- NormalizeData(cancer_MCF10, normalization.method = "LogNormalize", scale.factor = 1e5)
EP40V_MCF10 <- FindMarkers(cancer_MCF10, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V_MCF10 <- EP40V_MCF10[ ,c(2,5)]
colnames(EP40V_MCF10) <- c("MCF_log", "MCF_padj")
EP40V_MCF10$Gene <- rownames(EP40V_MCF10)
cancer_PC3 <- subset(cancer.integrated, subset=celltype=="PC3")
cancer_PC3 <- NormalizeData(cancer_PC3,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
cancer_PC3 <- NormalizeData(cancer_PC3, normalization.method = "LogNormalize", scale.factor = 1e5)
EP40V_PC3 <- FindMarkers(cancer_PC3, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
EP40V_PC3$lab <- rownames(EP40V_PC3)
EP40V_PC3 <- EP40V_PC3[ ,c(2,5)]
colnames(EP40V_PC3) <- c("PC3_log", "PC3_padj")
EP40V_PC3$Gene <- rownames(EP40V_PC3)
EP40V_PCMCF <- inner_join(x = EP40V_MCF10, y= EP40V_PC3, by = "Gene")
EP40V_PCMCF$bool <- FALSE
EP40V_PCMCF[EP40V_PCMCF$PC3_padj < 0.001 & EP40V_PCMCF$MCF_padj < 0.001, ]$bool <- TRUE
ggplot(EP40V_PCMCF, aes(x=MCF_log, y=PC3_log, label=Gene)) + geom_point()+
  geom_point(data=subset(EP40V_PCMCF, bool==TRUE), aes(x=MCF_log, y=PC3_log, label=Gene, color = "red"))+
  theme_classic()+ylim(-1, 1)
#
#Volcano
#
Diff <- subset(subset(subset(cancer_merge123456,subset=celltype=="MCF10A"), subset=run=="fifth", invert=TRUE), subset=condition=="normal")
Diff <- subset(subset(cancer.integrated, subset=celltype=="MCF10A"), subset=run=="fifth", invert=TRUE)
Diff <- subset(subset(subset(subset(cancer_merge123456,subset=celltype=="PC3"), subset=run==c("fifth", "sixth"), invert=TRUE), subset=condition=="normal"), subset=NEP=="75V", invert=TRUE)
Diff <- subset(subset(cancer.integrated, subset=celltype=="PC3"), subset=run==c("fifth", "sixth"), invert=TRUE)
DimPlot(Diff, group.by = "run")+DimPlot(Diff, group.by = "NEP")+DimPlot(Diff, group.by = "condition")+DimPlot(Diff, group.by = "celltype")
EP40V <- FindMarkers(Diff, group.by="NEP", ident.1 = "40V", min.pct = 0.25, logfc.threshold = 0)
rm(Diff)
EP40V$lab <- rownames(EP40V)
EP40V$lab[EP40V$p_val_adj >= 0.05] <- NA
ggplot(data=EP40V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP40V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  geom_text_repel()+theme_minimal()+theme_classic()
ggplot(data=EP40V, aes(x=avg_log2FC, y=-log10(p_val_adj), label = lab)) +
  geom_point(color = "black") + geom_point(data=subset(EP40V, p_val_adj < 0.05),aes(x=avg_log2FC, y=-log10(p_val_adj), color = "red"))+
  theme_minimal()+theme_classic()
