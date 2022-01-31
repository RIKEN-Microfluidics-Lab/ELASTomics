# WGCNA package
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
# "matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"
# "GO.db", "preprocessCore", "impute"
#BiocManager::install("WGCNA")
library(WGCNA)

options(stringsAsFactors = FALSE)
datExpr <- data.frame(t(GetAssayData(object=AML[["RNA"]])))
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


traitData <- t(GetAssayData(object=AML[["ADT"]]))


datTraits <- data.frame(traitData[,c(2:7)])
dim(datTraits)

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitColors),
                    main = "Sample dendrogram and trait heatmap")

#
# search gene network modules
source(file.path(rdir,"hunter_WGCNA_network_auto.R"))
#
# trait vs gene network modules
source(file.path(rdir,"hunter_WGCNA_relateModsToExt.R"))
#
# visualize the hierarchy of modules
source(file.path(rdir,"hunter_WGCNA_network_visualize.R"))



