
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 1, min.features = 1)

cellids <- colnames(pbmc)


dim(pbmc)
#head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)
rm(pbmc.data)

