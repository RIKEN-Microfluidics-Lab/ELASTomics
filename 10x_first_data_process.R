# load functions for barcode decoding
#source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
#source(file.path(rdir,'preprocess/preprocess_whitelist.R'))

# preprocess the count data and load reference
encode_barcode=FALSE
source(file.path(rdir,'preprocess/preprocess_RNAseq_data.R'))
#
# download reference data from ensembl with biomaRt
gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
source(file.path(rdir,'util/hunter_biomart_ref.R'))
#hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")
filter="ensembl_gene_id"
#symbol="mgi_symbol"
symbol="hgnc_symbol"
filter="ensembl_gene_id"
if (symbol=="mgi_symbol"){
  ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,filter,symbol))
}else{
  hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))
}

missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
# save count data with 10x format
rm(missing_ref,adding_ref)
source(file.path(rdir, 'preprocess/preprocess_save_10x_format.R'))

