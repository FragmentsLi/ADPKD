load('adpkd_all.RData')
library(Seurat)
library(SingleR)
library(celldex)

# Annotation
Idents(adpkd_qc) <- 'RNA_snn_res.0.7'
type<-c('PT','Tcell','MNP','MNP','Tcell',#0,1,2,3,4
        'FIB','Bcell','Tcell','PT','FIB',#5,6,7,8,9
        'PT','MAST','EPI','FIB','Tcell',#10,11,12,13,14
        'ENDO','EPI','EPI','Bcell','MNP',#15,16,17,18,19
        'Tcell','MNP','DT','FIB','EPI',#20,21,22,23,24
        'MAST','Bcell')#25,26

names(type) <- levels(adpkd_qc) 
adpkd_qc <- RenameIdents(adpkd_qc, type)       #修改Idents
adpkd_qc$celltype <- Idents(adpkd_qc)

# Cell markers
# Creating Figure 1C
data.usage<-DotPlot(adpkd_qc,feature=c('CD3D','CD3E','TRAC',#T Cells
                           'CD79A','MS4A1','CD19',#B Cells
                           'LYZ','CD14','CD68',#MNP
                           'KIT','ENPP3',#Mast Cells,KIT=CD117,ENPP3=CD203c
                           'PDGFRB','ACTA2','DCN','COL1A1','THY1',#Fibroblasts
                           'CDH1','EPCAM','KRT23','KRT18',#Epitheial Cells
                           'PECAM1','VWF',#Endotheial Cells
                           'SLC13A3','DCXR','GPX3',#PT Cells
                           'UMOD','DEFB1'#Distal tubule Cells
                          )#Distal tubule Cells
        ,group.by='celltype')$data

data.anno<-data.frame(features.plot=unique(data.usage$features.plot),
                      label=c(rep('TCell',3),rep('BCell',3),
                             rep('MNP',3),rep('MAST',2),
                             rep('FIB',5),rep('EPI',4),rep('ENDO',2),
                             rep('PT',3),rep('DT',2)))

df.plot<-plyr::join(data.usage,data.anno)

f.plot$id<-factor(df.plot$id,levels=sort(levels(df.plot$id)))

p<-ggplot(df.plot,aes(x=features.plot,y=id,size=pct.exp,color=avg.exp.scaled))+
geom_point()+scale_size("Percent Expressed",range=c(0,5))+
scale_color_gradientn(colours=viridis::viridis(20),
                      guide=guide_colorbar(ticks.colour='black',frame.colour='black'),
                      name="Average\nExpression",breaks=c(-0.5,1.0,2.5))+
cowplot::theme_cowplot()+
ylab("")+xlab("")+theme_bw()+
#scale_y_continuous(breaks=levels(df.plot$id))+
facet_grid(~label,scales="free_x",space="free")+theme_classic()+
theme(
  axis.text.x=element_text(size=8,angle=90,hjust=0.95,vjust=0.2,color='black'),
  axis.text.y=element_text(size=8,color='black'),
  axis.ticks.y=element_blank(),
  axis.text.y.right=element_blank(),
  axis.ticks.x=element_blank(),
  legend.title=element_text(size=8),
  legend.text=element_text(size=8),
  strip.text.x=element_text(size=8,face='bold',color='black',vjust=0.5,margin=margin(b=3,t=3)),
  strip.background=element_rect(colour='black',fill='white',size=1)
  )

png(filename = "ADPKD_all_features.png",width = 210, height = 70, 
    units = "mm", res = 300)
p
dev.off()

# Creating Figure 1B
library(Seurat)
library(ggsci)
library(ggplot2)
library(fplot)
adpkd_qc@meta.data$orig.ident2<-unlist(lapply(as.character(adpkd_qc@meta.data$orig.ident),function(x){
  return(paste0('Sample ',which(names(table(adpkd_qc@meta.data$orig.ident))==x)))}))
adpkd_qc@meta.data$orig.ident2<-ordered(adpkd_qc@meta.data$orig.ident2,
                          levels=c("Sample 1","Sample 2","Sample 3","Sample 4","Sample 5",
                                   "Sample 6","Sample 7","Sample 8","Sample 9","Sample 10"))

p<-DimPlot(adpkd_qc, reduction = 'umap',group.by = 'celltype',
        label = FALSE, label.size = 8,pt.size = 0.01) + scale_color_npg()+ labs(title=paste('Cell Types')) +
theme(axis.text=element_text(size=8),axis.title=element_text(size=8),
     legend.title=element_text(size=8),legend.text=element_text(size=8),
     plot.title=element_text(size=8))+
theme(legend.position="right")

png(filename = "ADPKD_all_UMAP_0.7_type.png",width = 120, height = 90, 
    units = "mm", res = 300)
p
dev.off()


