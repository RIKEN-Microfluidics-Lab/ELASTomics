# list of genes are obtained from https://www.nature.com/articles/nbt.4231#Sec32
# Paired-cell sequencing enables spatial gene expression mapping of liver endothelial cells Nat Biotech 2018
#
landmark_genes_cv <- read.csv('/home/samba/pihome/2021/Kuchimaru/genes_cv.csv',header = F)
landmark_genes_pn <- read.csv('/home/samba/pihome/2021/Kuchimaru/genes_pn.csv',header = F)
#
#

vivo <- subset(x=pbmc, subset=dish =="E00")

landmark_genes_cv <- subset(landmark_genes_cv, (V1 %in% ms_ref$gene_short_name))
landmark_genes_cv$type <- "pericentral"
landmark_genes_pn <- subset(landmark_genes_pn, (V1 %in% ms_ref$gene_short_name))
landmark_genes_pn$type <- "periportal"
landmark_genes_combine <- rbind(landmark_genes_cv,landmark_genes_pn)

landmark_genes_exp <- vivo[["RNA"]]@data[landmark_genes_combine$V1,]
rownames(landmark_genes_combine) <- landmark_genes_combine$V1
landmark_genes_type <- data.frame(landmark_genes_combine$type)
rownames(landmark_genes_type) <- landmark_genes_combine$V1

library(pheatmap)

pheatmap(landmark_genes_exp,
         annotation_row = landmark_genes_type,
         cluster_cols = TRUE,cluster_rows = TRUE,scale="row")
#
# ref https://www.nature.com/articles/nature21065
#
gse84498 <-read.table("/home/samba/pihome/2021/Kuchimaru/GSE84498_umitab.txt", header=TRUE, sep = "\t",row.names=1)
gene_list <- rownames(gse84498)
gene_list <- sapply(strsplit(gene_list, ";"), "[", 1)
rownames(gse84498)<-gene_list
gse84498_sub <- subset(x=gse84498, gene_list %in% ms_ref$gene_short_name)


# Initialize the Seurat object with the raw (non-normalized data).
gse84498seurat <- CreateSeuratObject(counts = gse84498, project = "gse84498", min.cells = 100, min.features = 1000)
gse84498seurat <- NormalizeData(gse84498seurat, normalization.method = "LogNormalize", scale.factor = 1e5)
gse84498seurat <- FindVariableFeatures(gse84498seurat, selection.method = "vst", nfeatures = 200)


all.genes <- rownames(gse84498seurat)
gse84498seurat <- ScaleData(gse84498seurat, features = all.genes)
gse84498seurat <- RunPCA(gse84498seurat, npcs=10, features = VariableFeatures(object = gse84498seurat))
gse84498seurat <- FindNeighbors(gse84498seurat, dims = 1:7)
gse84498seurat <- FindClusters(gse84498seurat, resolution = 0.1)

cluster = as.numeric(Idents(gse84498seurat))
gse84498seurat <- RunUMAP(gse84498seurat, dims = 1:4)
p1 <- DimPlot(gse84498seurat, reduction = "pca")
p2 <- DimPlot(gse84498seurat, reduction = "umap")
p1+p2
gse84498_landmark_genes <- gse84498_norm[landmark_genes_combine$V1,]

gse84498_hepatocyte <- subset(x=gse84498seurat, idents = c("0","1") )
DimPlot(gse84498_hepatocyte, reduction = "umap")
gse84498_norm<- gse84498_hepatocyte[["RNA"]]@data[landmark_genes_combine$V1,]


pheatmap(gse84498_landmark_genes,
         annotation_row = landmark_genes_type,
         cluster_cols = TRUE,cluster_rows = TRUE,scale="row")

p1 <- VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"))
tenx <- FetchData(vivo,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(vivo,vars="nFeature_RNA")
median(tenx$nFeature_RNA)

p2 <- VlnPlot(gse84498_hepatocyte, features = c("nCount_RNA","nFeature_RNA"))
tenx <- FetchData(gse84498_hepatocyte,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(gse84498_hepatocyte,vars="nFeature_RNA")
median(tenx$nFeature_RNA)
p1+p2
