library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Seurat)
library(ggplot2)
library(ggrepel)
exp.matrix <- data.frame(t(data.frame(tig.combined.nep[["integrated"]]@data)))
var.gene <- VariableFeatures(tig.combined.nep)

#exp.matrix <- as.matrix(exp.matrix[,200:500])
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
response <- data.frame(t(data.frame(tig.combined.nep[["DTD"]]@data)))

#cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)


#visualize correlation landscape with aging marker genes
cors$bool <- FALSE
cors[cors$gene %in% aging_gene,]$bool<-TRUE
#cors[cors$gene %in% "IGFBP2",]$bool<-TRUE
cors[cors$gene %in% anti_aging_genes,]$bool<-TRUE
#cors[cors$gene %in% gene_list_name$symbol,]$bool<-TRUE

ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors,pval<1e-10),aes(x=rank,y=cors, color = "red"))+
  geom_text_repel(data=subset(cors,bool==TRUE),aes(rank,cors,label=gene))+theme_bw()

#
# heatmap of correlated genes
#
cors.sig<-subset(subset(cors,subset=pval<1e-10),subset=abs(cors)>0.05)

exp.matrix <- data.frame(tig.combined.nep[["integrated"]]@data)
dtd <-FetchData(object = tig.combined.nep, vars = c("dtd_FLD500"))
rownames(dtd)<-colnames(exp.matrix)
age<-data.frame(tig.combined.nep[["age"]])
rownames(age)<-colnames(exp.matrix)
dtd$age <-age$age
out<-pheatmap(exp.matrix[rownames(cors.sig),order(dtd$dtd_FLD500,decreasing = FALSE)],
              cluster_cols = FALSE, show_colnames = FALSE,
              annotation_col = dtd,
              clustering_distance_rows = "correlation")

clust_genes<- sort(cutree(out$tree_row, k=3))
clust1_genes<-names(clust_genes[clust_genes==1])
clust2_genes<-names(clust_genes[clust_genes==2])
clust3_genes<-names(clust_genes[clust_genes==3])
#clust4_genes<-names(clust_genes[clust_genes==4])

gene.list<-symbol2entrez(clust3_genes)
ego_result <- enrichGO(gene          = gene.list$gene_id, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1,
                       qvalueCutoff  = 0.1, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE)
View(ego_result@result)
# FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "IER2",group.by = "age")+
#   FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "GNG11",group.by = "age")+
#   FeatureScatter(tig.combined.nep,feature1 = "dtd_FLD500",feature2 = "PAICS",group.by = "age")

clust_genes <- data.frame(clust_genes)
clust_genes$cors <- cors[rownames(clust_genes),]$cors
clust_genes$genes <- rownames(clust_genes)
clust_genes$bool <- FALSE
clust_genes[clust_genes$genes %in% c("ATF3","FOS","ADM","CDKN1A","DTYMK","NCL","PCLAF","DLGAP5"),]$bool<-TRUE
ggplot(clust_genes,aes(x=factor(clust_genes,levels=c("3","1","2")),y=cors))+geom_violin()+geom_jitter()+
  geom_label_repel(data=subset(clust_genes,bool==TRUE),
                   aes(x=factor(clust_genes,levels=c("3","1","2")),y=cors,label=genes))


#visualize correlation landscape with aging marker genes
cors$bool <- 0
cors[cors$gene %in% clust1_genes,]$bool<-1
cors[cors$gene %in% clust2_genes,]$bool<-2
cors[cors$gene %in% clust3_genes,]$bool<-3
#cors[cors$gene %in% gene_list_name$symbol,]$bool<-TRUE

ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors,bool==1),aes(x=rank+3,y=cors, colour = "1"))+
  geom_point(data=subset(cors,bool==3),aes(x=rank-3,y=cors, colour = "3"))+
  geom_point(data=subset(cors,bool==2),aes(x=rank-3,y=cors, colour = "2"))+
  geom_label_repel(data=subset(cors,bool==1),aes(rank,cors,label=gene))

#
# heatmap and clustering
#
library(reshape2)

clust1.exp <- exp.matrix[rownames(exp.matrix) %in% clust3_genes,]
clust1.sexp<- data.frame(apply(clust1.exp, 1, scale))
rownames(clust1.sexp)<-rownames(dtd)
clust1.sexp$dtd <-dtd$dtd_FLD500
clust1.sexp$age <-dtd$age
#View(melt(clust1.sexp,id.vars = c("dtd","age")))
clust.sexp<-melt(clust1.sexp,id.vars = c("dtd","age"))
p3<-ggplot(clust.sexp,aes(x=dtd,y=value))+
  geom_point()+
  geom_smooth()
p3+p1+p2

clust2.exp <- exp.matrix[rownames(exp.matrix) %in% clust2_genes,]
clust3.exp <- exp.matrix[rownames(exp.matrix) %in% clust3_genes,]

library(gplots)
senescence <-aging_gene
cor_gene <- cors.sig$gene
sasp<-c("YWHAH","SMC3","RAN","PAICS","NCL","MYH9","MSN","HSPD1","HNRNPDL","HNRNPAB","HMGA1","FABP5","CORO1C",
        "TIMP1","TFPI2","STOM","QSOX1","MFGE8","MFAP4","KLF4","IGFBP4","HLA.B","DCN")
data <- list(senescence = aging_gene, elast =cors.sig$gene,sasp=sasp)
venn(data)
