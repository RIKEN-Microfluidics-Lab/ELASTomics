library(stringr)
library(Seurat)
library(ggplot2)
library(R.utils)
library(reshape2)
library(gridExtra)
library("org.Hs.eg.db")
rdir <- "/home/samba/public/shintaku/github/ELASTomics/"
datadir <-"/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_TIGbulk/"
wdir <-"/home/samba/sanger/shintaku/ELASTomics/20220719HiSeqX014_TIGbulk/"

# laod whitelist and check the batch effect
source(file.path(rdir,"util/whitelist_encode.R"))
barcode <- read.table(file.path("/home/samba/sanger/shintaku/github/hunter/cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))
# create a list of whitelist files
filelist_whitelist <- data.frame(list.files(datadir,pattern="whitelist.txt"))
colnames(filelist_whitelist) <- c("filename")
# read whitelist file
for (icnt in 1:nrow(filelist_whitelist)){
  whitelist <- read.table(file.path(datadir,filelist_whitelist[icnt,]), sep="\t")
  colnames(whitelist) <- c('represent','variants', 'total','counts')
  # create functions for encoding
  encoded <- whitelist.umi_tools.list(whitelist,barcode$V1)
  encoded <- encoded[order(encoded$first_index,decreasing = FALSE),]
  correct_encoded <- encoded[encoded$lv_total<1,]
  correct_encoded$batch <- str_sub(filelist_whitelist[icnt,],1,10)
  rownames(correct_encoded) <- paste(correct_encoded$batch,correct_encoded$first_index,sep="-")
  if (icnt>1){
    allencoded <- rbind(allencoded,correct_encoded)
  }
  else{
    allencoded <- correct_encoded
  }
}
allencoded$GC <- as.numeric(lapply(lapply(as.character(allencoded$first_barcode),s2c),GC))
p0 <- ggplot(allencoded,aes(x=factor(GC),y=count))+ geom_violin()+ scale_y_log10()+geom_boxplot()
p1 <- ggplot(allencoded, aes(x = first_index, y = count, color=first_barcode))+geom_violin()
p2 <- ggplot(allencoded, aes(x = batch, y = count, color=batch))+geom_violin()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
grid.arrange(p0,p1, p2, nrow = 3)
rm(correct_encoded,encoded,p0,p1,p2)
active_barcode <- barcode[sort(unique(allencoded$first_index)),]
encode_barcode=TRUE

source(file.path(rdir,'unused/preprocess/preprocess_RNAseq_data.R'))

gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
source(file.path(rdir,'util/get_biomart_ref.R'))
filter="ensembl_gene_id"
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
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
# save count data with 10x format
rm(missing_ref,adding_ref)
source(file.path(rdir, 'unused/preprocess/preprocess_save_10x_format.R'))

tig.bulk <- Read10X(data.dir = wdir)

annotate.samples <-function(barcodes,PDL){
  cell <- barcodes
  cell$celltype <- PDL
  #  cell$culture <- culture
  colnames(cell)<-c("RTid","PDL")
  return(cell)
}

tig.pdl.33<-annotate.samples(data.frame(c("1","2","3")),"PDL33")
tig.pdl.44<-annotate.samples(data.frame(c("4","5","6")),"PDL44")
tig.pdl.52<-annotate.samples(data.frame(c("7","8","9")),"PDL52")
tig.pdl.57<-annotate.samples(data.frame(c("10","11","12")),"PDL57")
tig.pdl.65<-annotate.samples(data.frame(c("13","14","15")),"PDL65")
sample.meta.data <- rbind(tig.pdl.33,tig.pdl.44,tig.pdl.52,tig.pdl.57,tig.pdl.65)
rm(tig.pdl.33, tig.pdl.44, tig.pdl.52, tig.pdl.57, tig.pdl.65)
rownames(sample.meta.data)<-sample.meta.data$RTid

cellids <- colnames(tig.bulk) #substr(colnames(tig.bulk),12,22)

tig.bulk.seurat <- CreateSeuratObject(tig.bulk)
tig.bulk.seurat[["batch"]]<-substr(cellids,9,10)
tig.bulk.seurat[["RTid"]]<-substr(cellids,12,13)
tig.bulk.seurat[["PDL"]]<-sample.meta.data[substr(colnames(tig.bulk.seurat),12,13),]$PDL

FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "PDL")+
  FeatureScatter(tig.bulk.seurat,feature1="nCount_RNA",feature2 = "nFeature_RNA",group.by = "batch")

tig.bulk.seurat <- NormalizeData(tig.bulk.seurat, normalization.method = "LogNormalize", scale.factor = 1e5)
tig.bulk.seurat <- FindVariableFeatures(tig.bulk.seurat, selection.method = "vst", nfeatures = 1000)
all.genes <- rownames(tig.bulk.seurat)

tig.bulk.seurat <- ScaleData(tig.bulk.seurat, features = all.genes)
tig.bulk.seurat <- RunPCA(tig.bulk.seurat, npcs=3, features = VariableFeatures(object = tig.bulk.seurat))
DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "PDL")+
  DimPlot(tig.bulk.seurat, reduction = "pca",group.by = "batch")

tig.bulk.seurat <- subset(subset(tig.bulk.seurat,subset=batch=="03"), subset=PDL=="PDL57", invert= TRUE)
Idents(tig.bulk.seurat)<-tig.bulk.seurat[["PDL"]]

p1<-VlnPlot(tig.bulk.seurat, features = "RASD1",group.by = "PDL")+ylim(c(0,5))
p2<-VlnPlot(tig.bulk.seurat, features = "KLF2",group.by = "PDL")+ylim(c(0,5))
p3<-VlnPlot(tig.bulk.seurat, features = "RRAD",group.by = "PDL")+ylim(c(0,5))
p4<-VlnPlot(tig.bulk.seurat, features = "MFGE8",group.by = "PDL")+ylim(c(0,5))
p6<-VlnPlot(tig.bulk.seurat, features = "KLF4",group.by = "PDL")+ylim(c(0,5))
p7<-VlnPlot(tig.bulk.seurat, features = "FOS",group.by = "PDL")+ylim(c(0,5))
p8<-VlnPlot(tig.bulk.seurat, features = "BTG2",group.by = "PDL")+ylim(c(0,5))
p9<-VlnPlot(tig.bulk.seurat, features = "RGS3",group.by = "PDL")+ylim(c(0,5))
p5<-VlnPlot(tig.bulk.seurat, features = "SERTAD1",group.by = "PDL")+ylim(c(0,5))
p10<-VlnPlot(tig.bulk.seurat, features = "LINC02029",group.by = "PDL")+ylim(c(0,6))
gridExtra::grid.arrange(p1, p2,p3, p4, p5, p6,p7,p8,p9,p10,nrow = 2)
rm(p1, p2,p3, p4, p5, p6,p7,p8,p9,p10)

VlnPlot(tig.bulk.seurat, features = c("RRAD") ,group.by = "PDL", pt.size=0)+
  scale_fill_manual(values = c("#D55E00","#eab676","#1e81b0","#0072B2"))+NoLegend()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, , fill="black")

VlnPlot(tig.bulk.seurat, features = c("KLF2") ,group.by = "PDL", pt.size=0)+
  scale_fill_manual(values = c("#D55E00","#eab676","#1e81b0","#0072B2"))+NoLegend()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, , fill="black")

VlnPlot(tig.bulk.seurat, features = c("AC") ,group.by = "PDL", pt.size=0)+
  scale_fill_manual(values = c("#D55E00","#eab676","#1e81b0","#0072B2"))+NoLegend()+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, , fill="black")

FACS <- tig.bulk.seurat
FACS1 <-cbind(data.frame(FACS[["RNA"]]@data["RPL22L1",]), data.frame(FACS[["PDL"]]))
colnames(FACS1) <- c("RNA", "PDL")
amod <- aov(RNA ~ PDL, data = FACS1)
TukeyHSD(amod)



