# https://glmnet.stanford.edu/articles/glmnet.html
#https://www.r-bloggers.com/2021/05/lasso-regression-model-with-r-code/
library(glmnet)

exp.matrix <- t(data.frame(tig.combined.nep[["integrated"]]@data))
var.gene <- VariableFeatures(tig.combined.nep)
exp.matrix <- as.matrix(exp.matrix[,var.gene])
tig.combined.nep <- NormalizeData(tig.combined.nep,assay="DTD",normalization.method = "CLR",scale.factor = 1e2)
response <- data.frame(t(data.frame(tig.combined.nep[["DTD"]]@data)))

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


FeaturePlot(tig.combined.nep,features=c("dtd_FLD070","CDKN1A","FOS","JUN"))
vardata <-FetchData(tig.combined.nep,"FOS")
cor(vardata,response$FLD004, method="pearson")
  
