#
# hunter_SingleCellSignal
# http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html
# browseVignettes("SingleCellSignalR")


# Ligand/Receptor analysis using SingleCellSignalR
#signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)
# user guide
# https://rdrr.io/bioc/SingleCellSignalR/f/vignettes/UsersGuide.Rmd

# Visualization
#visualize(signal)
#intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)

library(SingleCellSignalR)

#scdata <- data.frame(GetAssayData(object=pbmc[["RNA"]]))
#scdata <- scdata-min(min(scdata))
# Data clustering
#all.genes <- rownames(pbmc)
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.1)


vitro1 <- subset(x=pbmc, subset=dish =="d01")
vitro2 <- subset(x=pbmc, subset=dish =="d02")
vitro3 <- subset(x=pbmc, subset=dish =="d03")
vitro4 <- subset(x=pbmc, subset=dish =="d04")
vitro <- merge(vitro1, y = vitro2, project = "vitro")
vitro <- merge(vitro, y = vitro3, project = "vitro")
vitro <- merge(vitro, y = vitro4, project = "vitro")
rm(vitro1,vitro2,vitro3,vitro4)


cluster = as.numeric(Idents(vitro))
scdata = data.frame(vitro[["RNA"]]@data)

all.genes <- row.names(scdata)
clust <- clustering(data=scdata, n.cluster=3, n=10,method="simlr",write=TRUE,pdf=FALSE)

signal = cell_signaling(data=scdata,genes=all.genes,cluster=clust$cluster,species ="mus musculus",
                        logFC=log2(1.0),s.score=0.5,int.type = "paracrine",write=TRUE,c.names=c("AMLaffected","E0771","AML"))

inter.net <- inter_network(data = scdata, signal = signal, genes = all.genes, cluster = clust$cluster, write = FALSE)
visualize_interactions(signal = signal)

vitro_result <- vitro[["dish"]]
vitro_result$cluster <- clust$cluster
vitro_result$tsne1 <- clust$`t-SNE`[,1]
vitro_result$tsne2 <- clust$`t-SNE`[,2]
vitro_result$normGFP <- vitro[["ADT"]]@data["Venus",]
vitro_result$HUNTER <- vitro[["RNA"]]@data["CS2nGFPT2AH2BmCherG01",]

ggplot(data=vitro_result,aes(x=tsne1,y=tsne2,color=dish))+geom_point()
