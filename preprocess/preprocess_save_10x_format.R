
# load reference data from biomaRt with filtering the gene_ids
# source('/home/samba/storage0/shintaku/HUNTER/hunter_biomart_ref.R')
allData$gene <- str_replace(str_replace(allData$gene,"mm10___",""),"GRCh38_","")
gene_list <- data.frame(str_replace(allData$gene,"_intron",""))
#ms_ref <- func.biomart.ref(ms_mart,gene_list,"mgi_symbol")
#ms_ref_intron <- func.ref.intron(ms_ref)
#all_ref <- rbind(ms_ref,ms_ref_intron)

# convert allData to gene matrix
gene_matrix <- tidyr::pivot_wider(allData, id_cols = gene, names_from = cell, values_from = count, values_fill=0)

# create barcodes.tsv.gz
barcodes <- data.frame(subset(colnames(gene_matrix),colnames(gene_matrix)!='gene'))
colnames(barcodes) <- "barcodes"
write.table(barcodes,file=paste0(wdir,"barcodes.tsv"), row.names=FALSE, col.names=FALSE, quote = FALSE)
gzip(paste0(wdir,'barcodes.tsv'),overwrite=TRUE)

# create features.tsv.gz
features <- data.frame(cbind(ms_ref$ensembl_gene_id,ms_ref$gene_short_name,"Gene Expression"))
colnames(features) <- c("ensembl_gene_id","gene_name","assay_type")
features <- data.frame(features[order(features$ensembl_gene_id),])
rownames(features) <- features$ensembl_gene_id
features_sel <- features[gene_matrix$gene,]

features_sel <- features_sel[which(!is.na(features_sel$ensembl_gene_id)),] # remove annotations absent in ref

write.table(features_sel,file=paste0(wdir,"features.tsv"), row.names=FALSE, sep = "\t", col.names=FALSE, quote = FALSE)
gzip(paste0(wdir,'features.tsv'),overwrite=TRUE)

# create matrix.mtx.gz
gene_matrix_mtx <- as.data.frame(gene_matrix[,2:ncol(gene_matrix)])
rownames(gene_matrix_mtx) <- gene_matrix$gene
gene_matrix_mtx <- gene_matrix_mtx[rownames(features_sel),]
gene_matrix_mtx <- Matrix(data.matrix(gene_matrix_mtx), sparse = TRUE)    # Thanks to Aaron for pointing this out
writeMM(gene_matrix_mtx,paste0(wdir,'matrix.mtx'))
gzip(paste0(wdir,'matrix.mtx'),overwrite=TRUE)

rm(allData,gene_list,gene_matrix,gene_matrix_mtx,barcodes,features,features_sel)

