library(ggplot2)
library(gridExtra)
library(IgorR)
library(Hmisc)
library(dplyr)
library(reshape2)
library(stringr)
library(minpack.lm)
library(exactRankTests)

#import file.ibw
datafile <- "/home/watson/pihome/2023/Shiomi/[AFM]/"     #file path
Day <- c("R050920_2", "R050920_3", "R051027", "R051107", "R051108", "R051108_2", "R051110", "R051110_2")  #measurement day
Lever <- c(5.75, 8.86, 6.10, 8.44, 8.44, 8.44, 8.25, 8.66)  #cantilever spring

slopeALL <- matrix(nrow = 140, ncol = 0)
slopeOriALL <- matrix(nrow = 2000, ncol = 0)
Young <- matrix(nrow = 0, ncol = 7)
colnames(Young) <- c("Wave","CellNo", "E","MT","R","base", "SD")
errorWave <- c()

for (d in 5:6){
  datadir <-paste0(datafile, Day[d], "/ForceCurve/")
  wdir <- paste0(datafile, Day[d], "/Result")
  CellTag <- read.csv(paste0(datafile, Day[d], "/CellTag.csv"), header = TRUE)
  Spring <- Lever[d]
  filelist_whitelist <- data.frame(list.files(datadir,pattern=".ibw"))
  colnames(filelist_whitelist) <- c("filename")
  #write.table(filelist_whitelist, file = "CellData.csv", append = F,sep = ",", row.names = F, quote = F) 
  for (icnt in 1:nrow(filelist_whitelist)){
    waveNo <- str_sub((filelist_whitelist[icnt,]), end = -5)
    CellNo <- CellTag[CellTag$Name==waveNo, 2]
    wavedata <- file.path(datadir,filelist_whitelist[icnt,])
    wavename <- read.ibw(wavedata, Verbose = FALSE, ReturnTimeSeries = FALSE, MakeWave = FALSE, HeaderOnly = FALSE)
    wavename <- as.data.frame(wavename)
    
    dt <- which.max(wavename[,1]) 
    colnames(wavename) <- c("move_um", "def_nm")
    wavename[,1] <- 1000000 * wavename[,1] 
    wavename[,2] <- 1000000000 * wavename[,2] 
    approach <- wavename[20:dt,1:2]
    approach[,1] <- approach[,1] - as.numeric(approach[1,1])
    approach[,2] <- approach[,2] - as.numeric(approach[1,2])
    #restract <- wavename[dt:nrow(wavename),1:2]
    rm(dt, wavedata, wavename)
    
    
    tryCatch({
      #Calculate baseline
      base <- round(nrow(approach) / 3)
      slopebase <- approach[10:base,]
      x <- slopebase[,1]
      y <- slopebase[,2]  
      fm <- nls(y~ s * x + t , start=c(s=0, t=0),)  
      fm_s <- coefficients(fm)["s"]
      fm_t <- coefficients(fm)["t"]
      approach$base_nm <- (fm_s * approach[,1] + fm_t)
      approach$line_nm <- as.matrix(approach["def_nm"] - approach["base_nm"])
      approach$sd_nm <- approach$line_nm * approach$line_nm
      
      #Calculate basepoint    
      fit <- c()
      for (t in 1:nrow(approach)) {
        if (t<20) {
          fit <- c(fit, 0)
        }else{
          fit <- c(fit, mean(approach$sd_nm[(t-19):t]))
        }
      }    
      approach$sd_nm <- fit    
      dev_p <-max(as.numeric(rownames(approach[(approach$sd_nm < (mean(approach$sd_nm[0:base])+ 2 * sd(approach$sd_nm[0:base]))), ])))  
      #dev_pa <-max(as.numeric(rownames(approach[(approach$line_nm < (mean(approach$line_nm[0:base])+ 4 * sd(approach$line_nm[0:base]))), ])))
      #dev_pb <-max(as.numeric(rownames(approach[(approach$line_nm < (mean(approach$line_nm[0:base])+ 2 * sd(approach$line_nm[0:base]))), ])))
      #dev_p <-dev_pb-(dev_pa - dev_pb) 
      dev <- approach[dev_p,1]                                                          #dev = start point
      dev_d <- approach$def_nm[dev_p]
      dev_t <- sum(approach$sd_nm[1:dev_p])/dev_p                                       #Difference from baseline
      curve <- approach[(dev < approach$move_um & approach$move_um < (dev + 1)),]
      curve <- curve[,c(-2, -3, -5)]
      
      
      #Calculate curcecurve
      curve$line_nm <- (curve$line_nm * Spring)                                         #cantilever spring
      colnames(curve) <- c("move_um", "force_pN")
      curve[,1] <- curve[,1] - as.numeric(curve[1,1])
      curve[,2] <- curve[,2] - as.numeric(curve[1,2])
      # #Sneddon model
      # pal <- (2/pi) * (1/(1-(0.5)^2)) * tan(pi/4)
      # x <- curve[,1]
      # y <- curve[,2]
      # fm <- nls(y~ pal * E * x^2, start=c(E=0), trace=TRUE)
      # fm_R <- modelr::rsquare(fm, curve)
      # fm_E <- coefficients(fm)["E"]
      # fm_MT <- fm_E
      # curve$Snedden <- (pal * fm_E * curve[,1] * curve[,1])
      # fit <- c()
      # for (t in 1:nrow(approach)) {
      #   if (t<dev_p) {
      #     fit <- c(fit, 0)
      #   }else{
      #     dis <- approach[t,1]-approach[dev_p,1]
      #     par <- as.numeric(pal * fm_E * dis * dis/ Spring)
      #     fit <- c(fit, par)
      #   }
      # }
      # approach$fitting <- fit
      
      
      #Yue Ding et al 2018
      curve_ed <- curve[-1,]
      x <- curve_ed$move_um
      obs <- curve_ed$force_pN
      curve <- curve[(curve$force_pN < 500),]
      tip_ang <- pi/4
      pred <- function(parS, xx) (2/pi) * tan(tip_ang) * parS$a * (x^2)*(1 + 0.95*(parS$b/(x*tan(tip_ang))))^0.92
      resid <- function(p, observed, xx) observed - pred(p,xx)
      parStart <- list(a = 1000, b = 1)
      nls.out <- nls.lm(par=parStart, fn=resid, observed=obs, xx=x, control=nls.lm.control(maxiter=1024,nprint=1))
      nls.out
      parEnd <- nls.out$par
      fm_E <- (1-(0.5^2)) * parEnd$a
      fm_s <- parEnd$b
      fm_MT <- parEnd$b * parEnd$a /2
      fm_R <- cor(obs, pred(nls.out$par, x))
      y <- pred(nls.out$par, x)
      fitted_data <- data.frame(x,y)
      
      
      #Data organization
      #g1 <- ggplot(approach, aes(x = move_um, y = sd_nm), colour="blue")+ geom_line()+
      #  geom_line(aes(move_um, base_nm), colour="red") +labs(title="Original")+
      #  geom_point(aes(dev, dev_d), colour="red") 
      #g2 <- ggplot() + geom_point(data = curve, mapping = aes(move_um, force_pN)) + 
      #  geom_line(data = fitted_data, mapping = aes(x,y)) + ggtitle(waveNo)+
      #  annotate(geom = "text",x=-Inf,y=Inf,hjust=-.2,vjust=2, label=paste("E = ", round(fm_E), "pN/um2",'\n', "MT = ", round(fm_MT), "pN/um"))
      #g <- gridExtra::grid.arrange(g1, g2, nrow = 2)
      #    dl_or_not <- askYesNo(msg='Do you want to save this data?', default = FALSE, prompts = 'y/n/na')
      #    if(dl_or_not==TRUE) {
      #      ggsave(paste("Res_", waveNo, ".pdf", sep=""), device = "pdf", width = 9, height = 9, path = wdir, g)
      Young <- rbind(Young, c(waveNo, CellNo, fm_E, fm_MT, fm_R, fm_s, dev_t))
      slopeALL <- cbind(slopeALL, curve[1:140,2])
      colnames(slopeALL)[ncol(slopeALL)] <- waveNo 
      slopeOriALL <- cbind(slopeOriALL, approach[1:2000,4])
      colnames(slopeOriALL)[ncol(slopeOriALL)] <- waveNo 
      #      }
    }, 
    error = function(e) {
      waveNo
      errorWave <- c(errorWave, waveNo)
    },
    finally = {
      message(icnt)
    },
    silent = TRUE
    )
  }
}

Young <- as.data.frame(Young)
Young[, 3] <- as.numeric(Young[, 3])
Young[, 4] <- as.numeric(Young[, 4])
Young[, 5] <- as.numeric(Young[, 5])
Young[, 6] <- as.numeric(Young[, 6])
Young[, 7] <- as.numeric(Young[, 7])
Setting <- c("Results/SD2_500nm_")
Young$Cell <- str_sub(Young$CellNo, start = 1, end = 5)
Young$CellID <- str_sub(Young$CellNo, start = 7)
Young$point <- str_sub(Young$Wave, start = -3, end = -2)
Res <- subset(Young, R>0.8 & E>1 & MT>1 & SD<10 & base<1)
Res <- group_by(Res, CellNo)
Res1 <- summarise(Res, n = n(), mean_E = mean(E), mean_ST = mean(MT), Med_E = median(E), Med_ST = median(MT), sd_E = sd(E), sd_ST = sd(MT))
Res1$Cell <- str_sub(Res1$CellNo, start = 1, end = 5)
Res1$CellLot <- str_sub(Res1$CellNo, start = 1, end = 8)
Res1 <- subset(Res1, Cell != "Dish")
ggplot(Res1, aes(x = Cell, y =Med_ST)) + geom_violin() + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 10) +theme_classic()
ggplot(Res1, aes(x = CellLot, y =Med_ST)) + geom_violin() + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 10) +theme_classic()
ggplot(Res1, aes(x = Med_ST, y = Med_E, color = Cell))+ geom_point()  +theme_classic()

Res <- group_by(Res, CellNo, point)
Res1 <- summarise(Res, n = n(), mean_E = mean(E), mean_ST = mean(MT))
Res1 <- group_by(Res1, CellNo)
Res1 <- summarise(Res1, n = n(), mean_E = mean(mean_E), mean_ST = mean(mean_ST), Med_E = median(mean_E), Med_ST = median(mean_ST), sd_E = sd(mean_E), sd_ST = sd(mean_ST))
Res1$Cell <- str_sub(Res1$CellNo, start = 1, end = 5)
Res1$CellLot <- str_sub(Res1$CellNo, start = 1, end = 8)
Res1 <- subset(Res1, Cell != "Dish")
#write.table(Res1, file = str_c(datafile, Setting, "CellData.csv", sep = ""), append = F,sep = ",", row.names = F, quote = F)
ggplot(Res1, aes(x = Cell, y =Med_ST)) + geom_violin() + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 4) +theme_classic()
ggplot(Res1, aes(x = Med_ST, y = Med_E, color = Cell))+ geom_point()  +theme_classic()



subset(Res1, subset=CellLot==c(""))

t.test(Med_ST ~ Cell, data = Res1)
t.test(Med_E ~ Cell, data = Res1)
wilcox.test(Res1[Res1$Cell == "PDL35",]$Med_ST, Res1[Res1$Cell == "PDL50",]$Med_ST)
t.test(Res1[Res1$Cell == "PDL35",]$Med_ST, Res1[Res1$Cell == "PDL50",]$Med_ST)
cor.test(CellInt$NorInt, CellInt$Med_ST)
cor(CellInt$NorDen, CellInt$Med_ST, method="spearman")

ggplot(curve, aes(x = move_um, y = force_pN))+ geom_point()  +theme_classic()

