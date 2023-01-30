library("org.Hs.eg.db")
library(clusterProfiler)

symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  gene.list <- gene.list[!duplicated(gene.list$symbol),]
  #geneset <- genemap[!duplicated(genemap[,1]), 2]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[gene.list$symbol,7] #correlation
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
gene_list_log2fc_cor <- symbol2entrez.order(MCF7)  
gse_result_cor<- gseGO(geneList     = gene_list_log2fc_cor,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)

symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
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
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 12,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.5,
                       verbose      = FALSE)

View(gse_result_cor@result)
View(gse_result_fol@result)

dotplot(gse_result_cor,showCategory = 12)
dotplot(gse_result_fol,showCategory = 12)
GSEA_cor <- gse_result_cor@result
GSEA_cor <- GSEA_cor[GSEA_cor$p.adjust<0.01, ]
GSEA_cor <- GSEA_cor[,c(1, 2, 4)]
colnames(GSEA_cor) <- c("ID", "Description", "Cor_Score")
GSEA_fol <- gse_result_fol@result
GSEA_fol <- GSEA_fol[GSEA_fol$p.adjust<0.01, ]
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

#cytoplasmic translation
g1 <- gseaplot(gse_result_fol, by = "runningScore", geneSetID = "GO:0002181")
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0002181")
gridExtra::grid.arrange(g1, g2, nrow = 2) 

#actin filament depolymerization
g1 <- gseaplot(gse_result_fol, by = "runningScore", geneSetID = "GO:0030042")
g2 <- gseaplot(gse_result_cor, by = "runningScore", geneSetID = "GO:0030042")
gridExtra::grid.arrange(g1, g2, nrow = 2)
