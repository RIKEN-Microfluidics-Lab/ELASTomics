# https://glmnet.stanford.edu/articles/glmnet.html
#https://www.r-bloggers.com/2021/05/lasso-regression-model-with-r-code/
library(glmnet)

exp.matrix <- data.frame(t(data.frame(tig.combined.nep[["integrated"]]@data)))
var.gene <- data.frame(VariableFeatures(tig.combined.nep))
response <- data.frame(t(data.frame(tig.combined.nep[["DTD"]]@data)))
#response <- FetchData(tig.combined.nep,"RRAD",slot="data")
rownames(response) <- rownames(exp.matrix)
cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors_p_val <- data.frame(t(apply(exp.matrix,2,function(x){unlist(cor.test(response$FLD004, x, method="pearson"))})))
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
cors$pval <- as.numeric(cors_p_val$p.value)
cors <- na.omit(cors)
cors.sig<-subset(subset(cors,subset=pval<0.001),subset=abs(cors)>0.15)
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_text(aes(label=ifelse(cors>0.2,as.character(gene),'')),hjust=0, vjust=0)
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]

ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_point(data=subset(cors, abs(cors) > 0.15), color = "red")+NoLegend()+theme_bw()+
#  geom_text_repel(data=subset(cors,gene %in% aging_gene),fontface = "italic", hjust=0, vjust=0)+
  geom_text_repel(data=subset(cors,gene %in% c('DLGAP5', 'PCLAF', 'NCL', 'DTYMK')),fontface = "italic", hjust=0, vjust=0)+NoLegend()+theme_bw()

#
# elastic net
#
exp.matrix <- t(data.frame(tig.combined.nep[["integrated"]]@data))
exp.matrix <- exp.matrix[,colnames(exp.matrix) %in% cors.sig[,]$gene]
alpha <- seq(0, 1, 0.01)
mse.df <- NULL
for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],mse = min(m$cvm)))
}
best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = best.alpha)
best.lambda <- m$lambda.min
en.model <- glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian",
                   lambda = best.lambda, alpha = best.alpha)
en.model.beta <-data.frame(coef(en.model, s="lambda.min"))
en.model.beta$gene <- rownames(en.model.beta)
en.model.beta$cors <- cors.sig[en.model.beta$gene,]$cors
en.model.beta <- na.omit(en.model.beta)
ggplot(en.model.beta, aes(x = cors, y = s1, label = gene))+geom_point()+theme_bw() +
  geom_point(data=subset(en.model.beta, abs(s1) > 0.1 & abs(cors - 0.05) > 0.25), color = "red")+
  geom_text_repel(aes(label=ifelse(abs(s1) > 0.1 & abs(cors - 0.05) > 0.25,as.character(gene),'')))

en.model.plus.beta <- en.model.beta[en.model.beta$s1>0,]
en.model.minus.beta <- en.model.beta[en.model.beta$s1<0,]
  
ggplot(en.model.beta,aes(x=cors,y=s1))+geom_point()+
  geom_point(data=subset(en.model.beta,gene %in% c("AC007952.4","AC091271.1","ABL2","FOS","RRAD","RHOB","RPL22L1","YWHAH","RASD1","AL021155.5","BTG2","KLF2")),aes(y=s1,x=cors,color="red"))+
  geom_text_repel(data=subset(en.model.beta,gene %in% c("AC007952.4","AC091271.1","ABL2","FOS","RRAD","RHOB","RPL22L1","YWHAH","RASD1","AL021155.5","BTG2","KLF2")),aes(y=s1,x=cors,label=gene,fontface = "italic"))+
  theme_bw()+ylim(c(-0.25,0.5))+NoLegend()
