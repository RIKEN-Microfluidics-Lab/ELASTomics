# https://glmnet.stanford.edu/articles/glmnet.html
#
library(glmnet)
library(coefplot)
exp.matrix <- t(data.frame(tig.nep[["RNA"]]@data))
var.gene <- VariableFeatures(tig)
exp.matrix <- exp.matrix[,var.gene]
tig.nep <- NormalizeData(tig.nep,assay="DTD",normalization.method = "CLR",scale.factor = 1)
response <- data.frame(t(data.frame(tig.nep[["DTD"]]@data)))

#
# elastic net
#
alpha <- seq(0, 1, 0.1)
mse.df <- NULL

for (i in 1:length(alpha)) {
  m <- cv.glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian", alpha = alpha[i])
  mse.df <- rbind(mse.df, data.frame(alpha = alpha[i],
                                     mse = min(m$cvm)))
}

best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]

m <- cv.glmnet(x = exp.matrix, y = response$FLD500, family = "gaussian", alpha = best.alpha)

best.lambda <- m$lambda.min

en.model <- glmnet(x = exp.matrix, y = response$FLD004, family = "gaussian",
                   lambda = best.lambda, alpha = best.alpha)

en.model.beta <-data.frame(en.model$beta)
#write.csv2(rownames(en.model.nonzero.beta),"/home/samba/public/shintaku/ELASTomics/elast.net.genes.fld004.csv")
#DoHeatmap(tig.nep,features = rownames(en.model.nonzero.beta))
#FeaturePlot(tig.nep,features=c("dtd_FLD004","RASD1","RPL22L1"))
#FeatureScatter(tig.nep,feature1="dtd_FLD500",feature2="CFAP100")
