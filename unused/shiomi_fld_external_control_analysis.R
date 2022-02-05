
#add meta data
FLDmapALL$plate <- str_sub(rownames(FLDmapALL),1,3)
FLDmapALL$cell <- str_sub(rownames(FLDmapALL),4,5)
FLDmapALL$treat <- str_sub(rownames(FLDmapALL),6,6)
FLDmapALL$gate <- str_sub(rownames(FLDmapALL),7,8)

# visualize
FLDsub <- subset(FLDmapALL,toupper(rownames(FLDmapALL)) %in% cellids)

FLDmelt <- melt(FLDsub, measure.vars=c("FLD010","FLD070","FLD500"))
ggplot(FLDmelt,aes(x=FLD004,y=value,color=variable))+geom_point()+scale_x_log10()+scale_y_log10()

FLDmelt <- melt(FLDsub, measure.vars=c("FLD004","FLD010","FLD070","FLD500"))
ggplot(FLDmelt,aes(x=factor(GC),y=value))+geom_boxplot()+scale_y_log10()


FLDsub <- FLDmapALL[FLDmapALL$plate=="T01",]
FLDmelt <- melt(FLDsub, measure.vars=c("FLD004","FLD010","FLD070","FLD500")) 
ggplot(FLDmelt,aes(x=FLDcon,y=value,color=variable))+geom_point()+scale_x_log10()+scale_y_log10()
