fucci_cellcycle_genes <- function(hs_ref){
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

cell_cycle_markers <- dplyr::left_join(cell_cycle_genes,hs_ref,by=c("geneID"="ensembl_gene_id"))
cell_cycle_markers <- cell_cycle_markers[!is.na(cell_cycle_markers$gene_short_name),]

# cell_cycle_markers_macosko_a <- read.csv("/home/samba/storage0/shintaku/macosco_cellcycle_genes.csv")
# cell_cycle_markers_macosko <- dplyr::left_join(cell_cycle_markers_macosko_a,all_ref,by=c("CCNE2"="gene_short_name"))
# cell_cycle_markers_macosko <- cell_cycle_markers_macosko[!is.na(cell_cycle_markers_macosko$ensembl_gene_id),]

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>% 
  dplyr::filter(phase == "S") %>% 
  pull("gene_short_name")
s_genes <- s_genes[!is.na(s_genes)]
# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_short_name")
g2m_genes <- g2m_genes[!is.na(g2m_genes)]

return(list(cell_cycle_markers,s_genes,g2m_genes))
}



