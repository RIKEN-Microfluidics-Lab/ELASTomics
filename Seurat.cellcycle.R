
source(file.path(rdir,"elast.biomaRt.R"))
source(file.path(rdir,"util/fucci_cellcycle_genes.R"))

sub_ref <- ref %>%
  dplyr::filter(gene_short_name %in% rownames(tig))

genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]


tig.bulk.seurat<- CellCycleScoring(tig.bulk.seurat,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
DimPlot(tig.bulk.seurat,reduction="pca",group.by = "celltype")+
DimPlot(tig.bulk.seurat,reduction="pca")


tig$CC.Difference <- tig$S.Score - tig$G2M.Score
#tig<- ScaleData(tig, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))

tig<- ScaleData(tig, vars.to.regress = "CC.Difference", features = c(s_genes, g2m_genes))
#features=c("S.Score","G2M.Score")
#FeaturePlot(pbmc, features = features,reduction = "pca")
p1 <- FeatureScatter(tig,feature1 = "G2M.Score",feature2 = "dtd_FLD004")+scale_y_log10()
p2 <- FeatureScatter(tig,feature1 = "S.Score",feature2 = "dtd_FLD004")+scale_y_log10()
p1+p2

FeatureScatter(tig,feature1 = "S.Score",feature2 = "G2M.Score")

tig <- RunPCA(tig, npcs=50, features = VariableFeatures(object = tig))
DimPlot(tig,reduction="pca")
FeaturePlot(tig,features="dtd_FLD004",reduction="pca")
#pbmc<-tig

VlnPlot(tig,features='dtd_FLD004')+scale_y_log10(limits=c(10,1e5))

