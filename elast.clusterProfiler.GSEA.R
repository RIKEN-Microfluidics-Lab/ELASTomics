#
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
#
#install.packages("msigdbr")
#library(msigdbr)
#BiocManager::install("fgsea")
#
# gse
# 
symbol2entrez.order <- function(genes){
  SYBOL2EG<-as.data.frame(unlist(org.Hs.egSYMBOL2EG))
  gene.list <- SYBOL2EG[SYBOL2EG$symbol %in% rownames(genes),]
  gene.list<-gene.list[!is.na(gene.list$gene_id),]
  rownames(gene.list) <- gene.list$symbol
  gene.list$s0 <- genes[rownames(genes) %in% gene.list$symbol,]
  gene.list <- gene.list[,c(1,3)]
  gene_list_order<- unlist(gene.list$s0)
  names(gene_list_order) <-as.character(gene.list$gene_id)
  gene_list_order <- gene_list_order[order(gene.list$s0,decreasing = T)]
  return(gene_list_order)
}
genes <-en.model.nonzero.beta
gene_list_log2fc <- symbol2entrez.order(genes)
#
# try BP: biological process, CC: cellular component, or MF: molecular function
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "CC",
                   minGSSize    = 12,
                   pAdjustMethod = "BH",
                   verbose      = FALSE)

ridgeplot(gse_result,showCategory = 20)
#d <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
#compare_cluster_GO_emap <- enrichplot::pairwise_termsim(gse_result, semData = d,  method="Wang")
#emapplot(compare_cluster_GO_emap, showCategory = 8)

#
# kegg
#
kk <- gseKEGG(gene_list_log2fc, nPerm=10000)
ridgeplot(kk)
gseaplot2(kk, geneSetID =1, title = kk$Description[1])
kk <- gseMKEGG(gene_list_log2fc, nPerm=10000)
ridgeplot(kk)
#gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])
gseaplot2(kk, geneSetID =1, title = kk$Description[1])
#
#
#
#browseKEGG(kk,)
