library(monocle)
tig <- subset(x=pbmc,subset=treat=='TA',invert=TRUE)
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
#tig <- subset(x=tig, subset=cell=='T2E',invert=TRUE)
tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 500)

all.genes <- rownames(tig)
#tig$CC.Difference <- tig$S.Score - tig$G2M.Score
tig <- ScaleData(tig, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(tig))
tig <- RunPCA(tig, features = VariableFeatures(object = tig))

print(tig[["pca"]], dims = 1:2, nfeatures = 50)

# plot plates/dishes/gates/pools/rtid
p1 <- DimPlot(tig, reduction = "pca")
p2 <- DimPlot(tig, reduction = "pca",group.by = "treat")
p3 <- DimPlot(tig, reduction = "pca",group.by = "cell")
p3 <- FeaturePlot(tig, reduction = "pca", features = "FLD004")
p1+p2+p3

pbmc_monocle <- as.CellDataSet(tig)


pbmc_monocle <- estimateSizeFactors(pbmc_monocle)
pbmc_monocle <- estimateDispersions(pbmc_monocle)
disp_table <- dispersionTable(pbmc_monocle)
ordering_genes <- subset(disp_table, mean_expression >= 50)
colnames(ordering_genes) <- "gene_id"


pbmc_monocle <- setOrderingFilter(pbmc_monocle, ordering_genes)
pbmc_monocle <- reduceDimension(pbmc_monocle)

pbmc_monocle <- orderCells(pbmc_monocle)
p1<-plot_cell_trajectory(pbmc_monocle, color_by = "cell")
p2<-plot_cell_trajectory(pbmc_monocle, color_by = "gate")
p1+p2

#plot_genes_in_pseudotime(pbmc_monocle, color_by = "cell")

tig[["Pseudotime"]] <- pbmc_monocle$Pseudotime

#FeatureScatter(tig,feature1 = "Pseudotime",feature2 = "Venus",group.by = "cell")+scale_y_log10()
FeatureScatter(tig,feature1 = "Pseudotime",feature2 = "FLD004",group.by = "cell")+scale_y_log10()
VlnPlot(tig,features='FLDcon',group.by = "treat")

#FeatureScatter(pbmc,feature1 = "GC",feature2 = "FLD004",group.by = "cell")+scale_y_log10()
#FeatureScatter(pbmc,feature1="Venus",feature2 = "FLD010",group.by = "cell")+scale_y_log10()+scale_x_log10()

