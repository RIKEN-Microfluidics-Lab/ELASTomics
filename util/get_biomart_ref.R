library(biomaRt)

func.biomart.ref <- function(hs_mart, gene_list,filter,symbol){
  reference=getBM(attributes=c("ensembl_gene_id","description",symbol,"gene_biotype","chromosome_name","entrezgene_id"),
                  filters=filter,values=gene_list,mart=hs_mart)
  
  colnames(reference) <- c("ensembl_gene_id","description","gene_short_name","gene_biotype","chromosome_name","entrezgene_id")
  reference <- reference[!duplicated(reference$ensembl_gene_id),]
  row.names(reference) <- reference$ensembl_gene_id
  empty_logical <- which(reference$gene_short_name=="")
  reference$gene_short_name <- as.character(reference$gene_short_name)
  reference[empty_logical,]$gene_short_name <- as.character(reference[empty_logical,]$ensembl_gene_id)
  
  return(reference)
}
func.ref.intron <- function(hs_ref){
  hs_ref_intron <- hs_ref
  hs_ref_intron$ensembl_gene_id <- paste0(hs_ref$ensembl_gene_id,"_intron")
  rownames(hs_ref_intron)<- hs_ref_intron$ensembl_gene_id
  return(hs_ref_intron)
}

