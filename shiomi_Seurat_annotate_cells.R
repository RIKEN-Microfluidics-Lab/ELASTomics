#pbmc[['plates']] <- substr(cellids,1,3)
#pbmc[['treat']] <- substr(cellids,4,5)
#pbmc[['cell']] <- substr(cellids,4,6)
#pbmc[['gate']] <- substr(cellids,7,8)
#pbmc[['pool']] <- substr(cellids,10,10)
#pbmc[['rtid']] <- substr(cellids,12,13)
#
# create reference table with gene_short_name
source(file.path(rdir,'/util/hunter_biomart_ref.R'))
gene_list <- data.frame(rownames(pbmc))
colnames(gene_list) <- "gene"
# #hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")
filter="hgnc_symbol"
hs_mart <- useMart(biomart="ensembl", host="useast.ensembl.org",dataset="hsapiens_gene_ensembl")
symbol="hgnc_symbol"
ms_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))
missing_ref <- subset(gene_list,!(gene %in% ms_ref$gene_short_name))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(hs_mart,symbol,missing_ref,adding_ref,gene_list)
