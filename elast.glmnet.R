# https://glmnet.stanford.edu/articles/glmnet.html
#https://www.r-bloggers.com/2021/05/lasso-regression-model-with-r-code/
library(glmnet)

exp.matrix <- t(data.frame(tig.combined.nep[["integrated"]]@data))
var.gene <- VariableFeatures(tig.combined.nep)

exp.matrix <- as.matrix(exp.matrix[,200:500])
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
response <- data.frame(t(data.frame(tig.combined.nep[["DTD"]]@data)))

#cor(exp.matrix$FOSB, exp.matrix$FOS)
cors <-apply(exp.matrix,2,function(x){cor(response$FLD004,x)})
cors<-data.frame(cors)
cors$gene <-rownames(cors)
cors$rank <-rank(cors$cors)
ggplot(cors,aes(x=rank,y=cors,label=gene))+geom_point()+
  geom_text(aes(label=ifelse(cors>0.2,as.character(gene),'')),hjust=0, vjust=0)




#
# elastic net
#
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

en.model.plus.beta <- en.model.beta[en.model.beta$s1>0,]
en.model.minus.beta <- en.model.beta[en.model.beta$s1<0,]

tig.combined.nep <- RenameIdents(tig.combined.nep, `0` = "TIG1-50", `1` = "TIG1-20")
DimPlot(tig.combined, group.by = "condition")+DimPlot(tig.combined.nep)

FeaturePlot(tig.combined.nep,features=c("dtd_FLD004","dtd_FLD070","dtd_FLD500",
                                        "ATP2B1-AS1","FOS",
                                        "FOSB","DUSP6","DLGAP5"),reduction = "umap")+
  DimPlot(tig.combined.nep)
  
vardata <-FetchData(tig.combined.nep,"NRG1")
cor(vardata,response$FLD004, method="pearson")
  
