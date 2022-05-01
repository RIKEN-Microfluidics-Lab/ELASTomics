#
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
#
#install.packages("msigdbr")
#library(msigdbr)
#BiocManager::install("fgsea")
#
# gse
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
  #gene.list$s0 <- genes[gene.list$symbol,1]
  gene.list$s0 <- genes[gene.list$symbol,1]
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}

gene_list_log2fc <- symbol2entrez.order(cors)
#
# try BP: biological process, CC: cellular component, or MF: molecular function
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 12,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
ridgeplot(gse_result,showCategory = 20)
View(gse_result@result)
gseaplot(gse_result, geneSetID = 1, by = "runningScore", title = kk$Description[1])

gene_list <- unlist(strsplit(gse_result@result[gse_result@result$Description=="aging",]$core_enrichment,"/"))
SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
gene_list_name <- SYBOL2EG[SYBOL2EG$gene_id %in% gene_list,]
gene_list_name
#
# kegg
#
kk <- gseKEGG(gene_list_log2fc,
              organism = "hsa", keyType = "kegg",
              exponent = 1, minGSSize = 10,
              maxGSSize = 500, pvalueCutoff = 1,
              pAdjustMethod = "BH", verbose = TRUE,
              use_internal_data = FALSE, seed = FALSE)


ridgeplot(kk)
View(kk@result)
gseaplot(kk, geneSetID =11, title = kk$Description[11])
gene_list <- unlist(strsplit(kk@result[kk@result$Description=="NF-kappa B signaling pathway",]$core_enrichment,"/"))
SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
gene_list_name <- SYBOL2EG[SYBOL2EG$gene_id %in% gene_list,]
gene_list_name

kk <- gseMKEGG(gene_list_log2fc, nPerm=10000)
ridgeplot(kk)
#gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])
gseaplot(kk, geneSetID =1, title = kk$Description[1])
#
#
#
browseKEGG(kk,pathID="hsa04064")
library(pathview)
pathview(gene.data = gene_list_log2fc, pathway.id = "hsa04110",
         limit = list(gene=1, cpd=1))
