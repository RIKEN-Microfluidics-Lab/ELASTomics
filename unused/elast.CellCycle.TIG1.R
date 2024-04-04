
source(file.path(rdir,"elast.biomaRt.R"))
source(file.path(rdir,"util/fucci_cellcycle_genes.R"))


sub_ref <- ref %>%
  dplyr::filter(gene_short_name %in% rownames(tig.combined.nep))

genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
tig.combined.nep <-CellCycleScoring(tig.combined.nep,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)

DimPlot(tig.combined.nep,reduction="umap")+
  FeatureScatter(tig.combined.nep,feature1="S.Score",feature2 = "dtd_FLD004")+
  FeatureScatter(tig.combined.nep,feature1="G2M.Score",feature2 = "dtd_FLD004")
VlnPlot(tig.combined.nep,features="dtd_FLD040")
RidgePlot(tig.combined.nep,features = c("dtd_FLD040","dtd_FLD500"))

