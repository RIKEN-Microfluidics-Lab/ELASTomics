library(dplyr)
library(ggplot2)
library(Seurat)
library(stringr)

### Violinplot and Boxplot ###
#Set the location of the downloaded file
rdir <- "/home/samba/pihome/2022/Shiomi/FACS_data/"

#select the analysis data
datafile <- str_c(rdir, "TIG-1_PDL/", sep = "")           #Figure 3j
datafile <- str_c(rdir, "TIG-1_siRRAD/", sep = "")        #Figure 3k
datafile <- str_c(rdir, "TIG-1_2-DG/", sep = "")          #Figure 3l
datafile <- str_c(rdir, "HeLa_pulse/", sep = "")          #Supplementary Figure 3a
datafile <- str_c(rdir, "HeLa_voltage/", sep = "")        #Supplementary Figure 3b
datafile <- str_c(rdir, "PC-3_FITC-BSA/", sep = "")       #Supplementary Figure 3h right
datafile <- str_c(rdir, "MDA-MB-231_FITC-BSA/", sep = "") #Supplementary Figure 3i right
datafile <- str_c(rdir, "MCF7_FITC-BSA/", sep = "")       #Supplementary Figure 3j right
datafile <- str_c(rdir, "MCF10A_FITC-BSA/", sep = "")     #Supplementary Figure 3k right
datafile <- str_c(rdir, "mHSPCs_FITC-BSA/", sep = "")     #Supplementary Figure 3l right
datafile <- str_c(rdir, "TIG-1_FITC-BSA/", sep = "")      #Supplementary Figure 3m right
datafile <- str_c(rdir, "MDA-MB-231/", sep = "")          #Supplementary Figure 3n
datafile <- str_c(rdir, "OVCAR-3/", sep = "")             #Supplementary Figure 3o
datafile <- str_c(rdir, "CHO-K1/", sep = "")              #Supplementary Figure 3p
datafile <- str_c(rdir, "GEM-81/", sep = "")              #Supplementary Figure 3q
datafile <- str_c(rdir, "K562/", sep = "")                #Supplementary Figure 3r
datafile <- str_c(rdir, "HeLa_Drug/", sep = "")           #Supplementary Figure 5e
datafile <- str_c(rdir, "HeLa_Cholesterol/", sep = "")    #Supplementary Figure 5f
datafile <- str_c(rdir, "mHSPCs+TER119/", sep = "")       #Supplementary Figure 9q

filelist_whitelist <-data.frame(list.files(datafile,pattern=".csv"))
#csv name:[date]_[celltype or cellcondition]_[applied voltage]
for (icnt in 1:nrow(filelist_whitelist)){
  samp <- str_split(filelist_whitelist[icnt,], "_",  n = 3)
  FACSdir <-data.frame(t(read.csv(paste0(datafile, filelist_whitelist[icnt,]), header = TRUE)))
  pbmc <- CreateSeuratObject(counts = FACSdir)
  pbmc[["hashtag"]]<-samp[[1]][2]
  pbmc[["name"]]<-str_c(samp[[1]][2], "_", str_sub(samp[[1]][3], end = -5), sep = "")
  pbmc[["voltage"]]<-str_sub(samp[[1]][3], end = -5)
  pbmc[["tag"]]<-str_sub(samp[[1]][3], end = 3)
  if(icnt == 1) {
    FACS <- pbmc
  } else {
    FACS<-merge(FACS,y=pbmc)
  }
}
FACS <- subset(FACS, subset= BL1.A>100)
FACS[["RNA"]]@counts["BL1.A", ] <- log10(FACS[["RNA"]]@counts["BL1.A", ])
VlnPlot(FACS, features = c("BL1.A"), group.by = "name", pt.size = 0.0, slot = "count")+
  geom_boxplot(width = 0.1, color = "black", fill="white")  +NoLegend()

#two-tailed Student’s t-test
FACS1 <-cbind(data.frame(FACS[["RNA"]]@counts["BL1.A",]), data.frame(FACS[["name"]]))
colnames(FACS1) <- c("BL1", "Hash")
t.test(BL1 ~ Hash, data = FACS1)

#Tukey’s t-test
FACS1 <-cbind(data.frame(FACS[["RNA"]]@counts["BL1.A",]), data.frame(FACS[["name"]]))
colnames(FACS1) <- c("BL1", "Hash")
amod <- aov(BL1 ~ Hash, data = FACS1)
TukeyHSD(amod)



### Ridgeplot ###
#Set the location of the downloaded file
rdir <- "/home/samba/pihome/2022/Shiomi/FACS_data/"

#select the analysis data
datafile <- str_c(rdir, "PC-3_PI/", sep = "")           #Supplementary Figure 3h left
datafile <- str_c(rdir, "MDA-MB-231_PI/", sep = "")     #Supplementary Figure 3i left
datafile <- str_c(rdir, "MCF7_PI/", sep = "")           #Supplementary Figure 3j left
datafile <- str_c(rdir, "MCF10A_PI/", sep = "")         #Supplementary Figure 3k left
datafile <- str_c(rdir, "mHSPCs_PI/", sep = "")         #Supplementary Figure 3l left
datafile <- str_c(rdir, "TIG-1_PI/", sep = "")          #Supplementary Figure 3m left

filelist_whitelist <-data.frame(list.files(datafile,pattern=".csv"))
for (icnt in 1:nrow(filelist_whitelist)){
  samp <- str_split(filelist_whitelist[icnt,], "_",  n = 3)
  FACSdir <-data.frame(t(read.csv(paste0(datafile, filelist_whitelist[icnt,]), header = TRUE)))
  pbmc <- CreateSeuratObject(counts = FACSdir)
  pbmc[["name"]]<-str_c(samp[[1]][2], "_", str_sub(samp[[1]][3], end = 3), sep = "")
  if(icnt == 1) {
    FACS <- pbmc
  } else {
    FACS<-merge(FACS,y=pbmc)
  }
}
FACS <- subset(FACS, subset=BL3.A>1)
RidgePlot(FACS, features = c("BL3.A"), group.by = "name")+NoLegend()

#Calculate dead cell rate
sample <- "TIG-1PDL44_00V"  #select sample name
PI <- subset(FACS, subset=name==sample)
100* dim(subset(PI, subset= BL3.A > 750))[2] / dim(PI)[2]
