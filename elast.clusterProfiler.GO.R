
#
# GOenrichmentAnalysis (experimental)
#BiocManager::install("org.Mm.eg.db")
#library("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library(clusterProfiler)
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
#
symbol2entrez<-function(gene.symbol){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% gene.symbol,]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  return(gene.list)
}


gene.list<-symbol2entrez(en.model.minus.beta$gene)

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
#ego_result.simple<-simplify(ego_result)
head(as.data.frame(ego_result))
#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster

gene.plus.list<-symbol2entrez(en.model.plus.beta$gene)
gene.minus.list<-symbol2entrez(en.model.minus.beta$gene)
gene.list <- list(plus=gene.plus.list$gene_id,
                  minus=gene.minus.list$gene_id)
xx <- compareCluster(gene.list, fun="groupGO",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC")
#summary(xx)
dotplot(xx)
#clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
#clusterProfiler::cnetplot(ego_result, categoryS ize="pvalue")
goplot(ego_result.simple)

