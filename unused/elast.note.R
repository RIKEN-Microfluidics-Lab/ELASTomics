
# annotate cells
source(file.path(rdir,"shiomi_Seurat_annotate_cells.R"))

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'shiomi_Seurat_10x_technical.R'))

# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/public/Shiomi/ELASTomicsindex"
channel <- c("Events","FSC","SSC","Venus","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")
source(file.path(rdir,'/io/hunter_Seurat_load_adt_data.R'))
#
# load cite-seq-count=FLD data as pbmc.tag
#
source(file.path(rdir,'io/hunter_Seurat_load_fld_data.R'))
# annotate cells with HTO and create FLD total
source(file.path(rdir,'20210816HiSeqX004_annotate_condition.R'))
#
#pbmc[["FLD"]] <- CreateAssayObject(counts=new_pbmc.tag[c("FLD004","FLD010","FLD040","FLD070","FLD150","Exter1","Exter2","Exter3","Exter4","Unmapp"),])
# check FLD by visualizing results
source(file.path(rdir,'20210816HiSeqX004_visualize_condition.R'))
# add HTO to Seurat object
pbmc[["HTO"]] <- CreateAssayObject(counts=pbmc.tag[c("T20CTL","T50CTL","TAZCTL"),])
#pbmc <- NormalizeData(pbmc, normalization.method = "CLR", scale.factor = 1e5,assay = "HTO")

DimPlot(pbmc,reduction="pca")
#p1<-DimPlot(pbmc,reduction="umap")
FeaturePlot(pbmc,features=c("FLD500","T20CTL","T50CTL","TAZCTL"),reduction = "pca")

#
source(file.path(rdir,"shiomi_fld_external_control_analysis.R"))
# first overview
# analyze data with PCA and UMAP
# find clusters and marker genes
#source(file.path(rdir,"shiomi_Seurat_technicalcheck.R"))

# check cell cycle dependence 
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
# find marker genes
source(file.path(rdir,"shiomi_Seurat_Marker.genes.R"))
# compute pseudotime and order cells along the gene expression
source(file.path(rdir,"shiomi_Seurat_monocle_pseudotime.R"))

#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))
