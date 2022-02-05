umi_div <- 4^3
umi_count<- 1:10:101

# expected number of non detected umi
umi_labels <- umi_div*(1-exp(-(umi_count/umi_div)))
# normalized expected number of detected umi
umi_dup <- 1-(umi_count-umi_labels)/umi_count
# normalized umi count by capacity
umi_count_fraction <- umi_count/umi_div
umi_data <- data.frame(t(rbind(umi_count_fraction,umi_dup)))
colnames(umi_data)<-c("read","umi")
umi_data$condition<-"predict"
ggplot(umi_data,aes(x=umi_count_fraction,y=umi_dup))+geom_point() +scale_x_log10()


#
# load experimental data
datadir <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana"
Read_count <- load.fld(datadir,"read_count",barcode)
UMI_count <- load.fld(datadir,"umi_count",barcode)
#
exp_data <- data.frame(Read_count$`FLD004-AACGTGAT`)
exp_data$UMI_count <- UMI_count$`FLD004-AACGTGAT`
colnames(exp_data) <- c("read","umi")
#
#
# normalize the umi counts by capacity
exp_data$norm_read <- exp_data$read/4^8
exp_data$norm_umi <- exp_data$umi/exp_data$read
# remove na values
exp_data_filter <- exp_data[!is.na(exp_data$norm_umi),]
# 
umi_exp <- data.frame(t(rbind(exp_data_filter$norm_read,exp_data_filter$norm_umi)))
colnames(umi_exp) <- c("read","umi")
umi_exp$condition <- "experiment"
#
# create merged data for comparison
#
merge_umi <- rbind(umi_data,umi_exp)
ggplot(merge_umi,aes(x=read,umi,color=condition))+geom_point()+scale_x_log10(limits=c(1e-2,1))

rm(umi_div,umi_count,umi_labels,umi_dup,umi_count_fraction,Read_count,UMI_count)
rm(umi_data,exp_data)
                                                                             
