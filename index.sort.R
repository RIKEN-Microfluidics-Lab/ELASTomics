datadir <- "/home/samba/sanger/shintaku/ELASTomics/20220922HiSeqX015_index"
index <- Read10X(data.dir = datadir)
index.dtd <- Read10X(data.dir = file.path(datadir,"CITE-seq"))

index.dtd.subset <- index.dtd[,str_sub(colnames(index.dtd),9,21) %in% str_sub(colnames(index),9,21)]
