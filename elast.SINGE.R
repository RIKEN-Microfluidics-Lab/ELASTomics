library(Matrix)
library(R.matlab)
X<-data.frame(t(exp.matrix[order(dtd$dtd_FLD500,decreasing = FALSE),rownames(cors.sig)]))
ptime <-dtd[order(dtd$dtd_FLD500,decreasing = FALSE),]$dtd_FLD500
gene_list<-rownames(X)
dataout <-"/home/samba/public/shintaku/ELASTomics/"
writeMat(file.path(dataout,"elast_tig.mat"),X=as(as.matrix(X),"sparseMatrix"),ptime=t(ptime))
writeMat(file.path(dataout,"gene_list.mat"),gene_list=gene_list)

network<-read.table("/home/samba/public/shintaku/ELASTomics/singe.txt",
                header = TRUE,sep = ",")
library(igraph)
subnetwork <- subset(network, subset=Regulator=="AC245014.3" & SINGE_Score>4 )
gg<-graph_from_edgelist(t(rbind(subnetwork$Regulator,subnetwork$Target)),directed = T)
plot(gg,edge.width=subnetwork$SINGE_Score)

subnetwork <- subset(network, subset=Regulator=="RRAD" & SINGE_Score>0.005 )
gg<-graph_from_edgelist(t(rbind(subnetwork$Regulator,subnetwork$Target)),directed = T)

plot(gg,edge.width=subnetwork$SINGE_Score,
     vertex.color="white",
     vertex.size=20,
     vertex.label.cex=0.8, 
     vertex.label.family="Arial",vertex.label.font=3)
FeatureScatter(tig,feature1="RRAD",feature2="STOM")
