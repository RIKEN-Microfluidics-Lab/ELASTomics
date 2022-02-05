
tig <- JackStraw(tig, num.replicate = 100)
tig <- ScoreJackStraw(tig, dims = 1:20)
JackStrawPlot(tig, dims = 1:20)
ElbowPlot(tig)

tig <- FindNeighbors(tig, dims = 1:19)
tig <- FindClusters(tig, resolution = 0.3)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(tig))
tig <- RunUMAP(tig, dims = 1:19)
p1 <- DimPlot(tig, reduction = "pca",group.by = "plate")
p2 <- DimPlot(tig, reduction = "umap",group.by = "plate")
p3<-DimPlot(tig)
p1+p2+p3

#find marker genes in each cluster
tig.markers <- FindAllMarkers(tig, only.pos = FALSE )

p1 <- DimPlot(tig, reduction = "umap",group.by = "plate")
p2<-FeaturePlot(tig,features="Venus",slot="count")
p3<-DimPlot(tig)

cluster <-data.frame(as.numeric(Idents(tig)))
colnames(cluster) <- "cluster"
rownames(cluster) <- colnames(tig)
cluster$gate <- tig[['gate']]
cluster$plate <- tig[['plate']]


cluster_density <- cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::count(plate, name = 'count') %>%
  dplyr::group_by(plate) %>%
  mutate(freq=count/sum(count))


p4<-ggplot(cluster_density,aes(y=cluster,x=plate$plate,fill=count))+ 
  geom_tile()+ theme_classic()+
  scale_fill_gradientn(colours = c("white", "red"))+geom_text(aes(label = count)) 

p1+p2+p3+p4

rm(p1,p2,p3,plot1,plot2,count_summary,count_summary_mean,pca_topcells,top10,all.genes,tenx)