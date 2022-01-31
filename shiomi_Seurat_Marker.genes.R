library(pheatmap)
#
# 5Aza cells/RG cells
#
EP.markers <- FindMarkers(pbmc,logfc.threshold = 0.2,ident.1=colnames(subset(pbmc,subset=type=="E")))
perturbed_gene_EP <- subset(EP.markers, p_val_adj<0.001)
shared_genes_EP <- rownames(perturbed_gene_EP)
shared_genes_EP <- ms_ref[ms_ref$gene_short_name %in% shared_genes_EP,]
perturbed_expression <- pbmc[["RNA"]]@data[shared_genes_EP$gene_short_name,]
#perturbed_expression <- data.frame(t(perturbed_expression))
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("type","Amu")]],
         cluster_cols = FALSE,cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, scale='row')

LH.markers <- FindMarkers(pbmc,logfc.threshold = 0.2,ident.1=colnames(subset(pbmc,subset=Amu=="Hole")))
perturbed_gene_LH <- subset(LH.markers, p_val_adj<0.001)
shared_genes_LH <- rownames(perturbed_gene_LH)
shared_genes_LH <- ms_ref[ms_ref$gene_short_name %in% shared_genes_LH,]
perturbed_expression <- pbmc[["RNA"]]@data[shared_genes_LH$gene_short_name,]
#perturbed_expression <- data.frame(t(perturbed_expression))
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("type","Amu")]],
         cluster_cols = FALSE,cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, scale='row')

shared_genes<- intersect(rownames(perturbed_gene_EP),rownames(perturbed_gene_LH))
shared_genes <- ms_ref[ms_ref$gene_short_name %in% shared_genes,]
perturbed_expression <- pbmc[["RNA"]]@data[shared_genes$gene_short_name[251:316],]
perturbed_expression <- pbmc[["RNA"]]@data[c("CXCL1", "CXCL2", "CXCL3", "CXCL8", "TGFB2", "TNFAIP3", "FAM53C", "EGR2", "BMP2", "RRAD", "FOS", "DDIT3", "C11orf96", "ZNF451-AS1"),]

perturbed_expression <- as.data.frame(t(perturbed_expression))
perturbed_expression$Rad <- pbmctagR$pFLDALL
perturbed_expression <- perturbed_expression[order(perturbed_expression$Rad),]
perturbed_expression$Rad <- NULL
perturbed_expression <- as.data.frame(t(perturbed_expression))

pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("type","Amu")]], 
         cluster_cols = FALSE,cluster_rows = FALSE, show_colnames = FALSE, show_rownames = TRUE, scale='row')


FLDmarker_FLD <- tig[["FLD"]]@data["pFLDRatio",]
FLDmarker_FLD <- replace(FLDmarker_FLD, which(is.na(FLDmarker_FLD)), 0)
FLDmarker_genes<- matrix(nrow = 0, ncol = 2)
colnames(FLDmarker_genes) <-c("gene", "correlation")
icnt <- (length(rownames(tig[["RNA"]]@data)))
for (i in 1:icnt) {
  FLDmarker_RNA <- tig[["RNA"]]@data[i,]
  if (sum(FLDmarker_RNA) >0) {
    p <- c(rownames(pbmc[["RNA"]]@data)[i], cor(FLDmarker_FLD, FLDmarker_RNA))
    FLDmarker_genes<- rbind(FLDmarker_genes, p)
    rm(p, FLDmarker_RNA)
  }
}
max(FLDmarker_genes[,2])

FLDmarker_RNA <- tig[["RNA"]]@data["COL1A2",]
x <- data.frame(
  FLD = FLDmarker_FLD,
  RNA  = FLDmarker_RNA)

g <- ggplot(x, aes(x = FLD, y = RNA))
g <- g + geom_point()
plot(g)




perturbed_gene <- list(RG=rownames(perturbed_gene_RG),HEA=rownames(perturbed_gene_HEA))
p0 <- venn.diagram(perturbed_gene,filename =NULL, fill=c(2,3), alpha=0.4, lty=3)
grid::grid.newpage()
grid::grid.draw(p0)
#
# intersection of 5Aza and RG
#
shared_genes_RG_HEA<- intersect(rownames(perturbed_gene_HEA),rownames(perturbed_gene_RG))
shared_genes_RG_HEA <- ms_ref[ms_ref$gene_short_name %in% shared_genes_RG_HEA,]
perturbed_expression <- pbmc[["RNA"]]@data[shared_genes_RG_HEA$gene_short_name,]
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("gate","cell")]],
         cluster_cols = TRUE,cluster_rows = TRUE,scale='row')

#
# control cells
#
cntl_HeLa <- subset(x=pbmc,subset=cell!=c("HEA") )
elp.markers <- FindMarkers(cntl_HeLa,logfc.threshold = 0.1,ident.1=colnames(subset(cntl_HeLa,subset=gate=="RG")))
perturbed_gene_elp <- subset(elp.markers, p_val_adj<0.05)
perturbed_gene_elp
#

#
# intersection of 5Aza and RG HeLa
#
RG_marker_regardless_of_5Aza <-intersect(shared_genes_RG_HEA$gene_short_name,rownames(perturbed_gene_elp))
RG_marker_regardless_of_5Aza <-ms_ref[ms_ref$gene_short_name %in% RG_marker_regardless_of_5Aza,]
write.csv(RG_marker_regardless_of_5Aza,"/home/samba/pihome/2021/Shintaku/Shiomi/RG_marker_regardless_of_5Aza_out_of_15genes.csv")

perturbed_expression <- pbmc[["RNA"]]@data[RG_marker_regardless_of_5Aza$gene_short_name,]
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("gate","cell")]],
         cluster_cols = TRUE,cluster_rows = TRUE,scale='row')

FeatureScatter(subset(x=cntl_HeLa,subset=cell=="HEE"),feature1="S.Score",feature2 = "G2M.Score",group.by = "gate")

#
#
#
write.csv(perturbed_gene_HEA,file="/home/samba/pihome/2021/Shintaku/Shiomi/marker.genes_HEA.csv")
write.csv(perturbed_gene_RG,file="/home/samba/pihome/2021/Shintaku/Shiomi/marker.genes_RG.csv")
write.csv(shared_genes_RG_HEA,file="/home/samba/pihome/2021/Shintaku/Shiomi/shared_marker.genes_RG_HEA.csv")
