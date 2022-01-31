#
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
#
#install.packages("msigdbr")
#library(msigdbr)
#BiocManager::install("fgsea")
#
# gse
# 
gene_list <- function(perturbed_gene_HEA,ms_ref){
  perturbed_gene_HEA_ref <-match(rownames(perturbed_gene_HEA),ms_ref$gene_short_name)
  perturbed_gene_HEA$entrez <- ms_ref[perturbed_gene_HEA_ref,]$entrez_annotation
  perturbed_gene_HEA <- perturbed_gene_HEA[!is.na(perturbed_gene_HEA$entrez),]
  gene_list_log2fc<- unlist(perturbed_gene_HEA$avg_log2FC)
  names(gene_list_log2fc) <-as.character(perturbed_gene_HEA$entrez)
  gene_list_log2fc <- gene_list_log2fc[order(gene_list_log2fc,decreasing = T)]
  return(gene_list_log2fc)
}
gene_list_log2fc <- gene_list(elp.markers,ms_ref)
#
# try BP: biological process, CC: cellular component, or MF: molecular function
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "CC",
                   minGSSize    = 12,
                   pvalueCutoff = 0.4,
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
