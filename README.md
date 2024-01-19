# ELASTomics
ELASTomics (electroporation-based lipid bilayer assay for surface tension and transcriptomics) is an approach for profiling the physical properties of plasma membranes and gene expression of single cells. To use this package, you will need the R statistical computing environment and several packages available through Bioconductor and CRAN.


## ELASTomics for Cancer cell line
cancer.main.R outputs figures corresconding to Figs. 1b, c, e, g-l, and Supplementary Fig. 3c-g, 3s-u, 6a-c, 7b-i, 8a-d in the manuscript of Shiomi et al. (2024).
This program is linked to the following code:
- elast.load.elast.data.R
- cancer1_Seurat.clustering.R
- cancer2_Seurat.clustering.R
- cancer3_Seurat.clustering.R
- cancer4_Seurat.clustering.R
- cancer5_Seurat.clustering.R
- cancer6_Seurat.clustering.R
- cancer.integrated.R
- cancer.MCF10A.R
- cancer.GSEA.R


## ELASTomics for mHSPCs
hsc_main.R outputs figures corresconding to Figs. 2b-i, and Supplementary Fig. 9b-o, 9r-u, 10a-e in the manuscript of Shiomi et al. (2024).
This program is linked to the following code:
- elast.load.elast.data.R
- Seurat.clustering.R


## ELASTomics for TIG-1 cells
tig_main.R outputs figures corresconding to Figs. 3b-g, and Supplementary Fig. 12 in the manuscript of Shiomi et al. (2024).
This program is linked to the following code:
- elast.load.elast.data.R
- Seurat.clustering.R
- elast.integrated.R
- elast.glmnet.R
- cancer.GSEA.R


## Bulk RNA-sequencing
bulk.tig.senescence.level.R outputs figures corresconding to Figs. 3h,i, and Supplementary Fig. 11e-m in the manuscript of Shiomi et al. (2024).
This program is linked to the following code:
- util/whitelist_encode.R
- unused/preprocess/preprocess_RNAseq_data.R
- unused/preprocess/preprocess_save_10x_format.R


## AFM analysis
AFManalysis_TIG-1.R and AFManalysis_MCF10A.R outputs figures corresconding to Figs. 1d and Supplementary Fig. 4f-i, 11d in the manuscript of Shiomi et al. (2024).
