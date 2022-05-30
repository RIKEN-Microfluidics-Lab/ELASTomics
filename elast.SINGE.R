library(Matrix)
library(R.matlab)
X<-data.frame(t(exp.matrix[order(dtd$dtd_FLD500,decreasing = FALSE),rownames(cors.sig)]))
ptime <-dtd[order(dtd$dtd_FLD500,decreasing = FALSE),]$dtd_FLD
gene_list<-rownames(X)
dataout <-"/home/samba/public/shintaku/ELASTomics/"
writeMat(file.path(dataout,"elast_tig.mat"),X=as.matrix(X),ptime=t(ptime))
writeMat(file.path(dataout,"gene_list.mat"),gene_list=gene_list)


library(igraph)
